#pragma once

#include <memory>
#include <span>
#include <vector>

namespace racon {

class Sequence;
class Overlap;

class Data {
 public:
  Data() = delete;
  ~Data();

  Data(const Data& that) = delete;
  Data& operator=(const Data&) = delete;

  Data(Data&&) noexcept;
  Data& operator=(Data&&) noexcept;

  std::span<const std::unique_ptr<Sequence>> sequences() const;
  std::span<const std::unique_ptr<Sequence>> targets() const;
  std::span<const std::unique_ptr<Sequence>> queries() const;
  std::span<const std::unique_ptr<Overlap>> overlaps(
      std::uint32_t target_id) const;

 private:
  struct Impl;
  friend Data LoadData(const std::string& sequences_path,
                       const std::string& overlaps_path,
                       const std::string& targets_path, double error_threshold,
                       bool keep_all);

  explicit Data(std::unique_ptr<Impl>);

  std::unique_ptr<Impl> pimpl_;
};

Data LoadData(const std::string& sequences_path,
              const std::string& overlaps_path, const std::string& targets_path,
              double error_threshold, bool keep_all);

}  // namespace racon
