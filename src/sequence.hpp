#pragma once


#include <string_view>
#include <memory>
#include <string>
#include <vector>

namespace bioparser {
template <class T>
class FastaParser;

template <class T>
class FastqParser;
}  // namespace bioparser

namespace racon {

class Sequence;
std::unique_ptr<Sequence> CreateSequence(const std::string& name,
                                         const std::string& data);

class Sequence {
 public:
  ~Sequence() = default;

  std::string_view name() const { return name_; }

  std::string_view data() const { return data_; }

  std::string_view reverse_complement() const { return reverse_complement_; }

  std::string_view quality() const { return quality_; }

  std::string_view reverse_quality() const { return reverse_quality_; }

  void create_reverse_complement();

  void transmute(bool has_name, bool has_data, bool has_reverse_data);

  friend bioparser::FastaParser<Sequence>;
  friend bioparser::FastqParser<Sequence>;
  friend std::unique_ptr<Sequence> CreateSequence(const std::string& name,
                                                  const std::string& data);

 private:
  Sequence(const char* name, uint32_t name_length, const char* data,
           uint32_t data_length);
  Sequence(const char* name, uint32_t name_length, const char* data,
           uint32_t data_length, const char* quality, uint32_t quality_length);
  Sequence(const std::string& name, const std::string& data);
  Sequence(const Sequence&) = delete;
  const Sequence& operator=(const Sequence&) = delete;

  std::string name_;
  std::string data_;
  std::string reverse_complement_;
  std::string quality_;
  std::string reverse_quality_;
};

}  // namespace racon
