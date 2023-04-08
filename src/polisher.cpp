#include "polisher.hpp"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string_view>

#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "overlap.hpp"
#include "sequence.hpp"
#include "spoa/spoa.hpp"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "window.hpp"

namespace racon {

auto BindSegmentsToWindows(std::span<const std::unique_ptr<Sequence>> sequences,
                           std::span<Window> windows, const Overlap& ovlp) {
  const auto cigar = ovlp.cigar();
  const auto q_id = ovlp.q_id();

  const auto q_len = ovlp.q_length();

  bool found_first_match = false;
  int32_t q_curr =
      (ovlp.strand() ? (ovlp.q_length() - ovlp.q_end()) : ovlp.q_begin()) - 1;
  auto q_first = q_curr + 1;

  int32_t t_curr = ovlp.t_begin() - 1;
  auto t_first = t_curr;

  auto q_last = q_first, t_last = t_first;
  auto window_idx = static_cast<std::size_t>(std::distance(
      windows.begin(),
      std::upper_bound(windows.begin(), windows.end(),
                       static_cast<uint32_t>(t_curr + 1),
                       [](uint32_t val, const Window& win) -> bool {
                         return val < win.last();
                       })));

  auto add_interval = [&]() {
    const auto len = q_last - q_first;

    const auto data = !ovlp.strand()
                          ? sequences[q_id]->data().substr(q_first, len)
                          : sequences[q_id]->reverse_complement().substr(
                                q_len - q_last + 1, len);

    const auto quality = [&] {
      if (sequences[q_id]->reverse_quality().empty()) {
        return std::string_view{};
      }
      return !ovlp.strand() ? sequences[q_id]->quality().substr(q_first, len)
                            : sequences[q_id]->reverse_quality().substr(
                                  q_len - q_last + 1, len);
    }();

    windows[window_idx].AddLayer(data, quality, t_first, t_last);
  };

  for (uint32_t i = 0, j = 0; window_idx < windows.size() && i < cigar.size();
       ++i) {
    if (cigar[i] == 'M' || cigar[i] == '=' || cigar[i] == 'X') {
      uint32_t k = 0, num_bases = atoi(&cigar[j]);
      j = i + 1;
      while (k < num_bases && window_idx < windows.size()) {
        ++q_curr;
        ++t_curr;

        if (!found_first_match) {
          found_first_match = true;
          q_first = q_curr;
          t_first = t_curr;
        }
        q_last = q_curr;
        t_last = t_curr;
        if (t_last == static_cast<int32_t>(windows[window_idx].last()) - 1) {
          if (found_first_match) {
            add_interval();
          }
          found_first_match = false;
          ++window_idx;
        }

        ++k;
      }
    } else if (cigar[i] == 'I') {
      q_curr += atoi(&cigar[j]);
      j = i + 1;
    } else if (cigar[i] == 'D' || cigar[i] == 'N') {
      uint32_t k = 0, num_bases = atoi(&cigar[j]);
      j = i + 1;
      while (k < num_bases && window_idx < windows.size()) {
        ++t_curr;
        if (t_curr == static_cast<int32_t>(windows[window_idx].last()) - 1) {
          if (found_first_match) {
            add_interval();
          }
          found_first_match = false;
          ++window_idx;
        }
        ++k;
      }
    } else if (cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P') {
      j = i + 1;
    }
  }
}

struct Polisher::Impl {
  Impl() = delete;

  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

  Impl(Impl&&) = delete;
  Impl& operator=(Impl&&) = delete;

  Impl(PolisherConfig config, Data data)
      : config(config), data(std::move(data)) {}

  std::vector<std::unique_ptr<Sequence>> Polish() {
    auto dst = std::vector<std::unique_ptr<Sequence>>();
    dst.reserve(data.targets().size());

    auto function_timer = biosoup::Timer();

    function_timer.Start();
    const auto n_targets = data.targets().size();

    auto n_aligned = std::atomic_size_t(0);
    auto n_polished = std::atomic_size_t(0);

    auto report_ticket = std::atomic_size_t(0);
    auto report_state = [&function_timer, n_targets, &n_aligned, &n_polished,
                         &report_ticket]() -> void {
      auto const to_percent = [n_targets](double val) -> double {
        return 100. * (val / n_targets);
      };

      if (auto ticket = ++report_ticket; ticket == report_ticket) {
        fmt::print(
            stderr,
            "\r[camel::ErrorCorrect]({:12.3f}) aligned {:3.3f}% | polished "
            "{:3.3f}%",
            function_timer.Lap(), to_percent(n_aligned),
            to_percent(n_polished));
      }
    };

    function_timer.Start();
    tbb::parallel_for(
        size_t(0), data.targets().size(), [&](size_t target_idx) -> void {
          const auto& target_seq = data.targets()[target_idx];
          const auto target_data = target_seq->data();
          const auto target_qual = target_seq->quality();

          const auto target_overlaps = data.overlaps(target_idx);
          tbb::parallel_for(
              size_t(0), target_overlaps.size(), [&](size_t ovlp_idx) -> void {
                data.overlaps(target_idx)[ovlp_idx]->cigar(data.sequences());
              });
          ++n_aligned;
          report_state();

          std::vector<Window> windows;
          windows.reserve(target_data.length() / config.window_length + 1);
          for (size_t pos = 0; pos < target_data.length();) {
            const auto nxt =
                std::min(pos + config.window_length, target_data.length());

            const auto data = target_data.substr(pos, config.window_length);
            const auto qual =
                target_qual.empty()
                    ? std::string_view()
                    : target_qual.substr(pos, config.window_length);

            windows.emplace_back(pos, nxt, data, qual);
            pos = nxt;
          }

          for (const auto& ovlp_ptr : target_overlaps) {
            BindSegmentsToWindows(data.sequences(), windows, *ovlp_ptr);
          }

          std::atomic<size_t> n_polished_windows;
          std::vector<std::string> window_consensues(windows.size());
          tbb::parallel_for(
              size_t(0), windows.size(), [&](size_t window_idx) -> void {
                auto consensus_flag = windows[window_idx].GenerateConsensus(
                    GetAlignmentEngine(config.poa_cfg), config.trim);
                n_polished_windows += consensus_flag.second;
                window_consensues[window_idx] = std::move(consensus_flag.first);
              });

          const auto seq_len = std::transform_reduce(
              window_consensues.cbegin(), window_consensues.cend(), size_t(0),
              std::plus<size_t>(), std::mem_fn(&std::string::size));

          std::string consensus_seq;
          consensus_seq.reserve(seq_len);
          for (const auto& it : window_consensues) {
            consensus_seq += it;
          }

          const auto polished_ratio =
              static_cast<double>(n_polished_windows) / windows.size();

          if (config.include_unpolished || polished_ratio > 0) {
            std::string tags;
            tags += " LN:i:" + std::to_string(seq_len);
            tags += " RC:i:" + std::to_string(target_overlaps.size());
            tags += " XC:f:" + std::to_string(polished_ratio);

            dst.push_back(CreateSequence(std::string(target_seq->name()) + tags,
                                         consensus_seq));
          }

          ++n_polished;
          report_state();
        });

    return dst;
  }

  spoa::AlignmentEngine& GetAlignmentEngine(POAConfig const config) {
    auto& engine = alignment_engines.local();
    if (!engine) {
      engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW, config.match, config.mismatch, config.gap);
    }

    return *engine;
  }

  PolisherConfig config;
  Data data;

  tbb::enumerable_thread_specific<std::unique_ptr<spoa::AlignmentEngine>>
      alignment_engines;
};

Polisher::Polisher(Polisher&& that) noexcept : pimpl_(std::move(that.pimpl_)) {}

Polisher& Polisher::operator=(Polisher&& that) noexcept {
  pimpl_ = std::move(that.pimpl_);
  return *this;
}

Polisher::~Polisher() {}

Polisher::Polisher(PolisherConfig config, Data data)
    : pimpl_(std::make_unique<Impl>(config, std::move(data))) {}

std::vector<std::unique_ptr<Sequence>> Polisher::Polish() {
  return pimpl_->Polish();
}

}  // namespace racon
