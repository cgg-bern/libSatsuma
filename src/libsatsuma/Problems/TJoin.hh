//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Config/Export.hh>
#include <lemon/list_graph.h>

namespace Satsuma {

struct SATSUMA_EXPORT TJoin {
    using GraphT = lemon::ListGraph;
    using Node = typename GraphT::Node;
    using Edge = typename GraphT::Edge;
    template<typename T> using EdgeMap = typename GraphT::EdgeMap<T>;
    template<typename T> using NodeMap = typename GraphT::NodeMap<T>;

    using CostScalar = double;
    using CostMap = GraphT::EdgeMap<CostScalar>;
    using TMap = GraphT::NodeMap<bool>;
    using Solution = GraphT::EdgeMap<bool>;

    GraphT const &g;
    TMap const &t;
    CostMap const &cost;

    TJoin(GraphT const& _g, TMap const&_t, CostMap const &_cost)
        : g(_g), t(_t), cost(_cost)
    {}

};

struct SATSUMA_EXPORT TJoinResult {
    std::unique_ptr<TJoin::Solution> solution;
    TJoin::CostScalar cost;
};

} // namespace Satsuma
