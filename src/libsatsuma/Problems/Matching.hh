//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <lemon/list_graph.h>

namespace Satsuma {

/// Max-weighted perfect matching
struct Matching {
    using GraphT = lemon::ListGraph;
    using Node = typename GraphT::Node;
    template<typename T> using EdgeMap = typename GraphT::EdgeMap<T>;
    template<typename T> using NodeMap = typename GraphT::NodeMap<T>;
    using Edge = typename GraphT::Edge;
    using WeightScalar = int64_t;
    using CostScalar = WeightScalar;
    GraphT g;
    GraphT::EdgeMap<WeightScalar> weight {g};

    using Solution = GraphT::EdgeMap<bool>;

    WeightScalar cost(Solution const&) const;
    bool is_perfect(Solution const&sol) const;
};

struct MatchingResult {
    std::unique_ptr<Matching::Solution> solution;
    Matching::CostScalar weight;
};

} // namespace Satsuma

