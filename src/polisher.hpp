#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "dataset.hpp"

namespace racon {

struct POAConfig {
  int8_t match;
  int8_t mismatch;
  int8_t gap;
};

struct PolisherConfig {
  uint32_t window_length;
  double quality_threshold;
  bool trim;

  POAConfig poa_cfg;
};

class Sequence;
class Overlap;

class Polisher {
 public:
  Polisher() = delete;

  Polisher(const Polisher&) = delete;
  Polisher& operator=(const Polisher&) = delete;

  Polisher(Polisher&&) noexcept;
  Polisher& operator=(Polisher&&) noexcept;

  ~Polisher();

  Polisher(Dataset dataset, PolisherConfig config);
  std::vector<std::unique_ptr<Sequence>> Polish();

 private:
  struct Impl;
  std::unique_ptr<Impl> pimpl_;
};

}  // namespace racon
