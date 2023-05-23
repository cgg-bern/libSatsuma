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

    struct EdgeInfo {
        Node u, v;
        bool u_head, v_head;
        CostScalar cost = 0.;
        int lower = 0;
        int upper = inf();
    };
    inline Edge add_edge(EdgeInfo const &info)
    {
        auto e = BiFlowGraph::add_edge(info.u, info.v);
        u_head[e] = info.u_head;
        v_head[e] = info.v_head;
        cost[e] = info.cost;
        lower[e] = info.lower;
        upper[e] = info.upper;
        return e;
    }

};

struct BiMCFResult {
    std::unique_ptr<BiMCF::Solution> solution;
    BiMCF::CostScalar cost;
    BiMCF::FlowScalar max_flow; // maximal flow on an Bi-MCF edge // TODO: remove again? not sure if useful
};

} // namespace Satsuma
