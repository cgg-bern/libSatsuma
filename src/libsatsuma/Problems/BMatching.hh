//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <lemon/list_graph.h>
#include <limits>

namespace Satsuma {

/// Max-weighted perfect capacitated b-matching
struct BMatching
{
    using GraphT = lemon::ListGraph;
    using Node = typename GraphT::Node;
    using Edge = typename GraphT::Edge;
    template<typename T> using EdgeMap = typename GraphT::EdgeMap<T>;
    template<typename T> using NodeMap = typename GraphT::NodeMap<T>;
    using CostScalar = double;
    using DegreeScalar = int;

    GraphT g;
    GraphT::NodeMap<DegreeScalar> degree {g};
    GraphT::EdgeMap<CostScalar> weight {g};
    GraphT::EdgeMap<DegreeScalar> capacity {g};

    using Solution = GraphT::EdgeMap<DegreeScalar>;

    bool is_valid(Solution const&sol) const;

    static DegreeScalar inf() {return std::numeric_limits<DegreeScalar>::max();}
};

struct BMatchingResult {
    std::unique_ptr<BMatching::Solution> solution;
    BMatching::CostScalar weight;
};


} // namespace Satsuma
