#include "io.hpp"

#include <stdexcept>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "bioparser/mhap_parser.hpp"
#include "bioparser/paf_parser.hpp"
#include "bioparser/parser.hpp"
#include "bioparser/sam_parser.hpp"
#include "sequence.hpp"

namespace racon {

static bool isSuffix(std::string_view src, std::string_view suffix) {
  if (src.size() < suffix.size()) {
    return false;
  }
  return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::unique_ptr<bioparser::Parser<Sequence>> createSequenceParser(
    const std::string& sequences_path) {
  if (isSuffix(sequences_path, ".fasta") ||
      isSuffix(sequences_path, ".fasta.gz") ||
      isSuffix(sequences_path, ".fna") || isSuffix(sequences_path, ".fna.gz") ||
      isSuffix(sequences_path, ".fa") || isSuffix(sequences_path, ".fa.gz")) {
    return bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
        sequences_path);
  } else if (isSuffix(sequences_path, ".fastq") ||
             isSuffix(sequences_path, ".fastq.gz") ||
             isSuffix(sequences_path, ".fq") ||
             isSuffix(sequences_path, ".fq.gz")) {
    return bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(
        sequences_path);
  } else {
    throw std::runtime_error(
        "[racon] error: "
        "file " +
        sequences_path +
        " has unsupported format extension (valid extensions: "
        ".fasta, .fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
        ".fq, .fq.gz)!\n");
  }
}

std::unique_ptr<bioparser::Parser<Overlap>> createOverlapParser(
    const std::string& overlaps_path) {
  if (isSuffix(overlaps_path, ".mhap") || isSuffix(overlaps_path, ".mhap.gz")) {
    return bioparser::Parser<Overlap>::Create<bioparser::MhapParser>(
        overlaps_path);
  } else if (isSuffix(overlaps_path, ".paf") ||
             isSuffix(overlaps_path, ".paf.gz")) {
    return bioparser::Parser<Overlap>::Create<bioparser::PafParser>(
        overlaps_path);
  } else if (isSuffix(overlaps_path, ".sam") ||
             isSuffix(overlaps_path, ".sam.gz")) {
    return bioparser::Parser<Overlap>::Create<bioparser::SamParser>(
        overlaps_path);
  } else {
    throw std::runtime_error(
        "[racon::createPolisher] error: "
        "file " +
        overlaps_path +
        " has unsupported format extension (valid extensions: "
        ".mhap, .mhap.gz, .paf, .paf.gz, .sam, .sam.gz)!\n");
  }
}

}  // namespace racon
