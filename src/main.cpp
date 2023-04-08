#include <cstdlib>

#include "cxxopts.hpp"
#include "data.hpp"
#include "fmt/core.h"
#include "polisher.hpp"
#include "tbb/task_arena.h"
#include "version.h"

int main(int argc, char** argv) {
  auto options = cxxopts::Options(
      "racon2", "racon2 is a stand-alone read and assembly polishing tool");

  /* clang-format off */
  options.add_options("misc")
    ("t,threads", "number of threads",
      cxxopts::value<uint32_t>()->default_value("1"))
    ("e,error-threshold", 
     "maximum allowed error rate used for filtering overlaps",
      cxxopts::value<double>()->default_value("0.3"));
  options.add_options("flags")
      ("f,fragment",
       "fragment correction instead of conting polishing")
      ("no-trimming",
       "disable consensus trimming at window ends\n");
  options.add_options("info")
    ("h,help", "print help")
    ("v,version", "print version and quit");
  options.add_options("window arguments")
    ("w,window-length", "size of window on which POA is performed",
      cxxopts::value<uint32_t>()->default_value("200"))
    ("q,quality-threshold", 
     "threshold for average base quality of windows used in POA",
      cxxopts::value<double>()->default_value("10.0"))
    ("m,match", "score for matching bases",
      cxxopts::value<int8_t>()->default_value("3"))
    ("x,mismatch", "score for mismatching bases",
      cxxopts::value<int8_t>()->default_value("-5"))
    ("g,gap", "gap penalty (must be nagative)",
      cxxopts::value<int8_t>()->default_value("-4"));
  options.add_options("input")
    ("sequences", "query sequences in FASTA/FASTQ format",
      cxxopts::value<std::string>())
    ("overlaps", "overlap file in MHAP/PAF/SAM format",
      cxxopts::value<std::string>())
    ("targets", "target sequences in FASTA/FASTQ format",
      cxxopts::value<std::string>());
  /* clang-format on */

  try {
    options.parse_positional({"sequences", "overlaps", "targets"});
    const auto result = options.parse(argc, argv);

    if (result.count("help")) {
      fmt::print(stderr, "{}", options.help());
      return EXIT_SUCCESS;
    }
    if (result.count("version")) {
      /* clang-format off */
      fmt::print(stderr, "racon {}.{}.{}",
        racon_VERSION_MAJOR,
        racon_VERSION_MINOR, 
        racon_VERSION_PATCH
      );
      /* clang-format on */
      return EXIT_SUCCESS;
    }

    tbb::task_arena arena(result["threads"].as<uint32_t>());
    arena.execute([&] {
      auto polisher = racon::Polisher(
          /* polisher config */
          racon::PolisherConfig{
              .window_length = result["window-length"].as<uint32_t>(),
              .quality_threshold = result["quality-threshold"].as<double>(),
              .trim = result["no-trimming"].as<bool>(),

              .poa_cfg =
                  racon::POAConfig{.match = result["match"].as<int8_t>(),
                                   .mismatch = result["mismatch"].as<int8_t>(),
                                   .gap = result["gap"].as<int8_t>()}},

          /* polisher data */
          racon::LoadData(result["sequences"].as<std::string>(),
                          result["overlaps"].as<std::string>(),
                          result["targets"].as<std::string>(),
                          result["error-threshold"].as<double>(),
                          result["fragment"].as<bool>()));
    });

  } catch (const std::exception& e) {
    fmt::print(stderr, "{}", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
