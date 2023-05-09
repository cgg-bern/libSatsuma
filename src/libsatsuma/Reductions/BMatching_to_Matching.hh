//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BMatching.hh>
#include <libsatsuma/Problems/Matching.hh>

namespace Satsuma {

class BMatching_to_Matching
{
public:
    struct Config {
        std::unique_ptr<Matching::NodeMap<BMatching::Node>> *out_orig_node = nullptr;
        std::unique_ptr<Matching::NodeMap<int>> *out_node_num = nullptr;
        std::unique_ptr<Matching::NodeMap<BMatching::Edge>> *out_internode_edge = nullptr;
    };
    BMatching_to_Matching(BMatching const &bm, Config const& _config = default_config);
    Matching const& matching() const { return matching_;}
    BMatchingResult translate_solution(const MatchingResult &matching_result) const;
private:
    static const Config default_config;
    BMatching const &bmatching_;
    Matching matching_;
    Matching::EdgeMap<BMatching::Edge> orig_edge{matching_.g};
    const double costmul_ = 1LL << 20;
};

} // namespace Satsuma

