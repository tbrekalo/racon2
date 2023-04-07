#pragma once

#include <memory>
#include <span>
#include <vector>

namespace racon {

class Sequence;
class Overlap;

class Dataset {
 public:
  Dataset() = delete;
  ~Dataset();

  Dataset(const Dataset& that) = delete;
  Dataset& operator=(const Dataset&) = delete;

  Dataset(Dataset&&) noexcept;
  Dataset& operator=(Dataset&&) noexcept;

  std::span<const std::unique_ptr<Sequence>> sequences() const;
  std::span<const std::unique_ptr<Sequence>> targets() const;
  std::span<const std::unique_ptr<Sequence>> queries() const;
  std::span<const std::unique_ptr<Overlap>> overlaps(
      std::uint32_t target_id) const;

 private:
  struct Impl;
  friend Dataset LoadDataset(const std::string& sequences_path,
                             const std::string& overlaps_path,
                             const std::string& targets_path,
                             double error_threshold, bool keep_all);

  explicit Dataset(std::unique_ptr<Impl>);

  std::unique_ptr<Impl> pimpl_;
};

Dataset LoadDataset(const std::string& sequences_path,
                    const std::string& overlaps_path,
                    const std::string& targets_path, double error_threshold,
                    bool keep_all);

}  // namespace racon
