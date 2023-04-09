#include "data.hpp"

#include <algorithm>
#include <iterator>
#include <memory>

#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "io.hpp"
#include "overlap.hpp"
#include "sequence.hpp"

namespace racon {

static constexpr uint32_t kChunkSize = 1024 * 1024 * 1024;  // ~ 1GB

Data::Data(Data&& that) noexcept : pimpl_(std::move(that.pimpl_)) {}

Data::~Data() {}

Data& Data::operator=(Data&& that) noexcept {
  pimpl_ = std::move(that.pimpl_);
  return *this;
}

Data::Data(std::unique_ptr<Impl> pimpl) : pimpl_(std::move(pimpl)) {}

struct Data::Impl {
  std::vector<std::unique_ptr<Sequence>> sequences;
  std::vector<std::vector<std::unique_ptr<Overlap>>> target_overlaps;

  std::span<std::unique_ptr<Sequence>> targets_span;
  std::span<std::unique_ptr<Sequence>> queries_span;
};

std::span<const std::unique_ptr<Sequence>> Data::sequences() const {
  return pimpl_->sequences;
}
std::span<const std::unique_ptr<Sequence>> Data::targets() const {
  return pimpl_->targets_span;
}
std::span<const std::unique_ptr<Sequence>> Data::queries() const {
  return pimpl_->queries_span;
}

std::span<const std::unique_ptr<Overlap>> Data::overlaps(
    std::uint32_t target_id) const {
  return pimpl_->target_overlaps[target_id];
}

Data LoadData(const std::string& sequences_path,
              const std::string& overlaps_path, const std::string& targets_path,
              double error_threshold, bool keep_all) {
  biosoup::Timer timer;
  timer.Start();

  auto seqs_parser = racon::createSequenceParser(sequences_path);
  auto ovlp_parser = racon::createOverlapParser(overlaps_path);
  auto trgs_parser = racon::createSequenceParser(targets_path);

  std::vector<std::unique_ptr<Sequence>> sequences;
  std::vector<std::unique_ptr<Overlap>> all_overlaps;

  sequences = trgs_parser->Parse(-1);
  const auto n_targets = sequences.size();

  fmt::print(
      stderr,
      "[racon2::loadAndFormatDataset]({:12.3f}) loaded {} target sequences\n",
      timer.Stop(), n_targets);

  timer.Start();

  std::unordered_map<std::string, uint64_t> name_to_id;
  std::unordered_map<uint64_t, uint64_t> id_to_id;
  for (uint64_t i = 0; i < n_targets; ++i) {
    name_to_id[fmt::format("{}t", sequences[i]->name())] = i;
    id_to_id[i << 1 | 1] = i;
  }

  std::vector<char> has_name(n_targets, true);
  std::vector<char> has_data(n_targets, true);
  std::vector<char> has_reverse_data(n_targets, false);

  uint64_t sequences_size = 0;

  seqs_parser->Reset();
  while (true) {
    uint64_t l = sequences.size();
    auto reads = seqs_parser->Parse(kChunkSize);
    if (reads.empty()) {
      break;
    }

    sequences.insert(sequences.end(), std::make_move_iterator(reads.begin()),
                     std::make_move_iterator(reads.end()));

    uint64_t n = 0;
    for (uint64_t i = l; i < sequences.size(); ++i, ++sequences_size) {
      auto it = name_to_id.find(fmt::format("{}", sequences[i]->name()));
      if (it != name_to_id.end()) {
        if (sequences[i]->data().size() !=
                sequences[it->second]->data().size() ||
            sequences[i]->quality().size() !=
                sequences[it->second]->quality().size()) {
          throw std::runtime_error(
              fmt::format("[racon::Polisher::initialize] error: "
                          "duplicate sequence {} with unequal data\n",
                          sequences[i]->name()));
        }

        name_to_id[fmt::format("{}q", sequences[i]->name())] = it->second;
        id_to_id[sequences_size << 1 | 0] = it->second;

        sequences[i].reset();
        ++n;
      } else {
        name_to_id[fmt::format("{}q", sequences[i]->name())] = i - n;
        id_to_id[sequences_size << 1 | 0] = i - n;
      }
    }

    sequences.erase(
        std::stable_partition(std::next(sequences.begin(), l), sequences.end(),
                              [](const std::unique_ptr<Sequence>& seq) -> bool {
                                return static_cast<bool>(seq);
                              }),
        sequences.end());
  }

  if (sequences_size == 0) {
    throw std::runtime_error(
        "[racon::Polisher::initialize] error: "
        "empty sequences set!\n");
  }

  has_name.resize(sequences.size(), false);
  has_data.resize(sequences.size(), false);
  has_reverse_data.resize(sequences.size(), false);

  fmt::print(
      stderr,
      "[racon2::loadAndFormatDataset]({:12.3f}) loaded {} query sequences\n",
      timer.Stop(), sequences.size() - n_targets);

  timer.Start();
  auto remove_invalid_overlaps = [&](uint64_t begin, uint64_t end) -> void {
    for (uint64_t i = begin; i < end; ++i) {
      if (all_overlaps[i] == nullptr) {
        continue;
      }
      if (all_overlaps[i]->error() > error_threshold ||
          all_overlaps[i]->q_id() == all_overlaps[i]->t_id()) {
        all_overlaps[i].reset();
        continue;
      }
      if (!keep_all) {
        for (uint64_t j = i + 1; j < end; ++j) {
          if (all_overlaps[j] == nullptr) {
            continue;
          }
          if (all_overlaps[i]->length() >= all_overlaps[j]->length()) {
            all_overlaps[j].reset();
          } else {
            all_overlaps[i].reset();
            break;
          }
        }
      }
    }
  };

  ovlp_parser->Reset();
  uint64_t c = 0;
  while (true) {
    auto overlaps_chunk = ovlp_parser->Parse(kChunkSize);
    if (overlaps_chunk.empty()) {
      break;
    }
    all_overlaps.insert(all_overlaps.end(),
                        std::make_move_iterator(overlaps_chunk.begin()),
                        std::make_move_iterator(overlaps_chunk.end()));

    uint64_t l = c;
    for (uint64_t i = l; i < all_overlaps.size(); ++i) {
      all_overlaps[i]->transmute(sequences, name_to_id, id_to_id);

      if (!all_overlaps[i]->IsValid()) {
        all_overlaps[i].reset();
        continue;
      }

      while (all_overlaps[c] == nullptr) {
        ++c;
      }
      if (all_overlaps[c]->q_id() != all_overlaps[i]->q_id()) {
        remove_invalid_overlaps(c, i);
        c = i;
      }
    }

    uint64_t n = 0;
    for (uint64_t i = l; i < c; ++i) {
      if (all_overlaps[i] == nullptr) {
        ++n;
      }
    }
    c -= n;
    all_overlaps.erase(
        std::stable_partition(std::next(all_overlaps.begin(), l),
                              all_overlaps.end(),
                              [](const std::unique_ptr<Overlap>& ovlp) -> bool {
                                return static_cast<bool>(ovlp);
                              }),
        all_overlaps.end());
  }
  remove_invalid_overlaps(c, all_overlaps.size());
  all_overlaps.erase(std::stable_partition(
                         std::next(all_overlaps.begin(), c), all_overlaps.end(),
                         [](const std::unique_ptr<Overlap>& ovlp) -> bool {
                           return static_cast<bool>(ovlp);
                         }),
                     all_overlaps.end());

  for (const auto& it : all_overlaps) {
    if (it->strand()) {
      has_reverse_data[it->q_id()] = true;
    } else {
      has_data[it->q_id()] = true;
    }
  }

  std::unordered_map<std::string, uint64_t>().swap(name_to_id);
  std::unordered_map<uint64_t, uint64_t>().swap(id_to_id);

  if (all_overlaps.empty()) {
    throw std::runtime_error(
        "[racon::Polisher::initialize] error: "
        "empty overlap set!\n");
  }

  for (size_t j = 0; j < sequences.size(); ++j) {
    sequences[j]->transmute(has_name[j], has_data[j], has_reverse_data[j]);
  }

  fmt::print(stderr,
             "[racon2::loadAndFormatDataset]({:12.3f}) loaded {} overlaps\n",
             timer.Stop(), all_overlaps.size());

  std::span<std::unique_ptr<Sequence>> targets_span(
      sequences.begin(), std::next(sequences.begin(), n_targets));
  std::span<std::unique_ptr<Sequence>> queries_span(
      std::next(sequences.begin(), n_targets), sequences.end());

  std::vector<std::vector<std::unique_ptr<Overlap>>> overlaps(n_targets);
  for (size_t i = 0; i < all_overlaps.size(); ++i) {
    overlaps[all_overlaps[i]->t_id()].push_back(std::move(all_overlaps[i]));
  }

  return Data(std::make_unique<Data::Impl>(
      Data::Impl{.sequences = std::move(sequences),
                 .target_overlaps = std::move(overlaps),
                 .targets_span = targets_span,
                 .queries_span = queries_span}));
}
}  // namespace racon
