//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <lemon/list_graph.h>
#include <limits>

namespace Satsuma {

// Min-Cost-Flow
struct MCF
{
    using GraphT = lemon::ListDigraph;
    using Node = typename GraphT::Node;
    using Arc = typename GraphT::Arc;
    template<typename T> using ArcMap = typename GraphT::ArcMap<T>;
    template<typename T> using NodeMap = typename GraphT::NodeMap<T>;
    using CostScalar = int64_t;
    using FlowScalar = int;

    static FlowScalar inf() {return std::numeric_limits<FlowScalar>::max();}

    GraphT g;
    NodeMap<FlowScalar> supply {g}; // TODO: globally use either supply or demand! lemon uses supply.
    ArcMap<CostScalar> cost {g};
    ArcMap<FlowScalar> upper {g};
    ArcMap<FlowScalar> lower {g}; // TODO: WARNING: ignored for regular BiMDF pipeline and assumed to be 0!
    using Solution = ArcMap<FlowScalar>;

    CostScalar compute_cost(Solution const& sol) const;

    FlowScalar inflow(Node n, MCF::Solution const&sol) const;
    FlowScalar outflow(Node n, MCF::Solution const&sol) const;

    inline Node add_node(FlowScalar _supply) {
        auto nh = g.addNode();
        supply[nh] = _supply;
        return nh;
    }

    struct ArcInfo {
        Node u, v;
        CostScalar cost = 0;
        FlowScalar lower = 0;
        FlowScalar upper = inf();
    };

    inline Arc add_edge(ArcInfo const &info) {
        auto arc = g.addArc(info.u, info.v);
        cost[arc] = info.cost;
        upper[arc] = info.upper;
        lower[arc] = info.lower;
        return arc;
    }
};

struct MCFResult {
    std::unique_ptr<MCF::Solution> solution;
    MCF::CostScalar cost;
};


} // namespace Satsuma

