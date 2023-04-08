#pragma once

#include "bioparser/parser.hpp"
#include "overlap.hpp"
#include "sequence.hpp"

namespace racon {

std::unique_ptr<bioparser::Parser<Sequence>> createSequenceParser(
    const std::string& sequences_path);

std::unique_ptr<bioparser::Parser<Overlap>> createOverlapParser(
    const std::string& overlaps_path);

}  // namespace racon
