//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <lemon/list_graph.h>

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

    GraphT g;
    NodeMap<FlowScalar> supply {g}; // TODO: globally use either supply or demand! lemon uses supply.
    ArcMap<CostScalar> cost {g};
    ArcMap<FlowScalar> upper {g};
    ArcMap<FlowScalar> lower {g}; // TODO: WARNING: ignored for regular BiMDF pipeline and assumed to be 0!
    using Solution = ArcMap<FlowScalar>;

    CostScalar compute_cost(Solution const& sol) const;

    FlowScalar inflow(Node n, MCF::Solution const&sol) const;
    FlowScalar outflow(Node n, MCF::Solution const&sol) const;
};

struct MCFResult {
    std::unique_ptr<MCF::Solution> solution;
    MCF::CostScalar cost;
};


} // namespace Satsuma

