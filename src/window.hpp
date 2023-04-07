#pragma once

#include <memory>
#include <string_view>

namespace spoa {
class AlignmentEngine;
}

namespace racon {

class Window {
 public:
  Window(uint32_t first, uint32_t last, std::string_view backbone,
         std::string_view quality);

  ~Window();

  Window(Window&&) noexcept;
  Window& operator=(Window&&) noexcept;
  
  /*
   *  @brief second params is true if polish is success,
   *  false in case the coverage was to low
   * */
  std::pair<std::string, bool> GenerateConsensus(
      spoa::AlignmentEngine& alignment_engine, bool trim);

  void AddLayer(std::string_view sequence, std::string_view quality,
                uint32_t first, uint32_t last);

  std::uint32_t first() const { return first_; }
  std::uint32_t last() const { return last_; }

 private:
  struct Impl;

  std::uint32_t first_;
  std::uint32_t last_;
  std::unique_ptr<Impl> pimpl_;
};

}  // namespace racon
