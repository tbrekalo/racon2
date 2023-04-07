#pragma once

#include <memory>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace bioparser {
template <class T>
class MhapParser;

template <class T>
class PafParser;

template <class T>
class SamParser;
}  // namespace bioparser

namespace racon {

class Sequence;

class Overlap {
 public:
  ~Overlap();

  uint32_t q_id() const { return q_id_; }
  uint32_t q_begin() const { return q_begin_; }
  uint32_t q_end() const { return q_end_; }
  uint32_t q_length() const { return q_length_; }

  uint32_t t_id() const { return t_id_; }
  uint32_t t_begin() const { return t_begin_; }
  uint32_t t_end() const { return t_end_; }
  uint32_t t_length() const { return t_length_; }

  uint32_t strand() const { return strand_; }

  bool IsValid() const { return is_valid_; }

  void Transmute(const std::vector<std::unique_ptr<Sequence>>& sequences,
                 const std::unordered_map<std::string, uint64_t>& name_to_id,
                 const std::unordered_map<uint64_t, uint64_t>& id_to_id);

  uint32_t length() const { return length_; }
  double error() const { return error_; }

  std::string_view cigar(std::span<const std::unique_ptr<Sequence>> sequences);
  std::string_view cigar() const { return cigar_; };

 private:
  friend bioparser::MhapParser<Overlap>;
  friend bioparser::PafParser<Overlap>;
  friend bioparser::SamParser<Overlap>;

  // #ifdef CUDA_ENABLED
  //   friend class CUDABatchAligner;
  // #endif

  Overlap(uint64_t a_id, uint64_t b_id, double accuracy, uint32_t minmers,
          uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
          uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);
  Overlap(const char* q_name, uint32_t q_name_length, uint32_t q_length,
          uint32_t q_begin, uint32_t q_end, char orientation,
          const char* t_name, uint32_t t_name_length, uint32_t t_length,
          uint32_t t_begin, uint32_t t_end, uint32_t matching_bases,
          uint32_t overlap_length, uint32_t maping_quality);
  Overlap(const char* q_name, uint32_t q_name_length, uint32_t flag,
          const char* t_name, uint32_t t_name_length, uint32_t t_begin,
          uint32_t mapping_quality, const char* cigar, uint32_t cigar_length,
          const char* t_next_name, uint32_t t_next_name_length,
          uint32_t t_next_begin, uint32_t template_length, const char* sequence,
          uint32_t sequence_length, const char* quality,
          uint32_t quality_length);
  Overlap();
  Overlap(const Overlap&) = delete;
  const Overlap& operator=(const Overlap&) = delete;

  std::string q_name_;
  uint64_t q_id_;
  uint32_t q_begin_;
  uint32_t q_end_;
  uint32_t q_length_;

  std::string t_name_;
  uint64_t t_id_;
  uint32_t t_begin_;
  uint32_t t_end_;
  uint32_t t_length_;

  uint32_t strand_;
  uint32_t length_;
  double error_;
  std::string cigar_;

  bool is_valid_;
  bool is_transmuted_;
  std::vector<std::pair<uint32_t, uint32_t>> breaking_points_;
};

}  // namespace racon
