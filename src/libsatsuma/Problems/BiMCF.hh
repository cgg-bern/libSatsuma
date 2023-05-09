//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiFlowGraph.hh>

namespace Satsuma {

/// Minimum-cost flow in bidirected graphs
struct BiMCF : public BiFlowGraph
{
    using GraphT = BiFlowGraph::GraphT;
    using Node = BiFlowGraph::Node;
    using Edge = BiFlowGraph::Edge;
    using FlowScalar = BiFlowGraph::FlowScalar;

    template<typename T> using EdgeMap = typename BiFlowGraph::EdgeMap<T>;
    template<typename T> using NodeMap = typename BiFlowGraph::NodeMap<T>;

    using CostScalar = double;

    EdgeMap<CostScalar> cost {g};
    CostScalar compute_cost(Solution const& sol) const;
};

struct BiMCFResult {
    std::unique_ptr<BiMCF::Solution> solution;
    BiMCF::CostScalar cost;
    BiMCF::FlowScalar max_flow; // maximal flow on an Bi-MCF edge // TODO: remove again? not sure if useful
};

} // namespace Satsuma
