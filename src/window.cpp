#include "window.hpp"

#include <algorithm>
#include <cstdlib>
#include <utility>

#include "fmt/core.h"
#include "spoa/spoa.hpp"

namespace racon {

struct Window::Impl {
  std::pair<std::string, bool> GenerateConsensus(
      spoa::AlignmentEngine& alignment_engine, bool trim) {
    std::string consensus;
    if (sequences.size() < 3) {
      return {std::string(sequences.front().begin(), sequences.front().end()),
              false};
    }

    /* clang-format off */
    spoa::Graph graph{};
    if (qualities.front().empty()) {
      graph.AddAlignment(spoa::Alignment(),
          sequences.front().data(), sequences.front().length());
    } else {
      graph.AddAlignment(spoa::Alignment(),
          sequences.front().data(), sequences.front().length(),
          qualities.front().data(), qualities.front().length());
    }
    /* clang-format on */

    std::vector<uint32_t> rank;
    rank.reserve(sequences.size());
    for (uint32_t i = 0; i < sequences.size(); ++i) {
      rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
      return positions[lhs].first < positions[rhs].first;
    });

    uint32_t offset = 0.01 * sequences.front().length();
    for (uint32_t j = 1; j < sequences.size(); ++j) {
      uint32_t i = rank[j];

      spoa::Alignment alignment;
      if (positions[i].first < offset &&
          positions[i].second > sequences.front().length() - offset) {
        alignment = alignment_engine.Align(sequences[i].data(),
                                           sequences[i].length(), graph);
      } else {
        std::vector<const spoa::Graph::Node*> mapping;
        auto subgraph =
            graph.Subgraph(positions[i].first, positions[i].second, &mapping);
        alignment = alignment_engine.Align(sequences[i].data(),
                                           sequences[i].length(), subgraph);
        subgraph.UpdateAlignment(mapping, &alignment);
      }

      /* clang-format off */
      if (qualities[i].empty()) {
        graph.AddAlignment(
          alignment,
          sequences[i].data(), sequences[i].length());
      } else {
        graph.AddAlignment(
          alignment,
          sequences[i].data(), sequences[i].length(),
          qualities[i].data(), qualities[i].length());
      }
      /* clang-format on */
    }

    std::vector<uint32_t> coverages;
    consensus = graph.GenerateConsensus(&coverages);

    if (trim) {
      uint32_t average_coverage = (sequences.size() - 1) / 2;

      int32_t begin = 0, end = consensus.size() - 1;
      for (; begin < static_cast<int32_t>(consensus.size()); ++begin) {
        if (coverages[begin] >= average_coverage) {
          break;
        }
      }
      for (; end >= 0; --end) {
        if (coverages[end] >= average_coverage) {
          break;
        }
      }

      if (begin >= end) {
        // TODO: log stuff
      } else {
        consensus = consensus.substr(begin, end - begin + 1);
      }
    }

    return {consensus, true};
  }

  void AddLayer(std::string_view sequence, std::string_view quality,
                std::uint32_t first, std::uint32_t last) {
    if (sequence.empty() || first == last) {
      return;
    }

    if (!quality.empty() && sequence.length() != quality.length()) {
      throw std::runtime_error(
          "[racon::Window::AddLayer] error: "
          "unequal quality size!\n");
    }

    if (first >= last || first > sequences.front().length() ||
        last > sequences.front().length()) {
      throw std::runtime_error(
          fmt::format("[racon::Window::AddLayer] error: "
                      "layer begin and end positions are invalid!\n (first, "
                      "last) = ({}, {})",
                      first, last));
    }

    sequences.emplace_back(sequence);
    qualities.emplace_back(quality);
    positions.emplace_back(first, last);
  }

  std::vector<std::string_view> sequences;
  std::vector<std::string_view> qualities;
  std::vector<std::pair<uint32_t, uint32_t>> positions;
};

Window::Window(uint32_t begin_pos, uint32_t end_pos, std::string_view backbone,
               std::string_view quality)
    : first_(begin_pos), last_(end_pos), pimpl_(std::make_unique<Impl>()) {
  pimpl_->sequences.push_back(backbone);
  pimpl_->qualities.push_back(quality);
  pimpl_->positions.emplace_back(0, 0);
}

Window::~Window() {}

Window::Window(Window&& that) noexcept
    : first_(std::exchange(that.first_, 0)),
      last_(std::exchange(that.last_, 0)),
      pimpl_(std::move(that.pimpl_)) {}

Window& Window::operator=(Window&& that) noexcept {
  this->first_ = std::exchange(that.first_, 0);
  this->last_ = std::exchange(that.first_, 0);
  this->pimpl_ = std::move(that.pimpl_);

  return *this;
}

void Window::AddLayer(std::string_view sequence, std::string_view quality,
                      uint32_t begin, uint32_t end) {
  pimpl_->AddLayer(sequence, quality, begin - first_, end - first_);
}

std::pair<std::string, bool> Window::GenerateConsensus(
    spoa::AlignmentEngine& alignment_engine, bool trim) {
  return pimpl_->GenerateConsensus(alignment_engine, trim);
}

}  // namespace racon
