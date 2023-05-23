//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiFlowGraph.hh>
#include <libsatsuma/Problems/CostFunction.hh>
#include <memory>
#include <variant>
#include <cmath>

namespace Satsuma {


/// Minimum-deviation flow in bidirected graphs
struct BiMDF : public BiFlowGraph
{
    using GraphT = BiFlowGraph::GraphT;
    using Node = BiFlowGraph::Node;
    using Arc = BiFlowGraph::Arc;
    using Edge = BiFlowGraph::Edge;
    using FlowScalar = BiFlowGraph::FlowScalar;
    using Solution = BiFlowGraph::Solution;
    using Guess = BiFlowGraph::Guess;
    template<typename T> using EdgeMap = typename BiFlowGraph::EdgeMap<T>;
    template<typename T> using NodeMap = typename BiFlowGraph::NodeMap<T>;

    using TargetScalar = double;
    using CostScalar = double;

    EdgeMap<CostFunction::Function> cost_function {g};

    struct EdgeInfo {
        Node u, v;
        bool u_head, v_head;
        CostFunction::Function cost_function = CostFunction::Zero();
        FlowScalar lower = 0;
        FlowScalar upper = inf();
    };


    inline Edge add_edge(EdgeInfo const &info)
    {
        auto e = BiFlowGraph::add_edge(info.u, info.v);
        u_head[e] = info.u_head;
        v_head[e] = info.v_head;
        cost_function[e] = info.cost_function;
        lower[e] = info.lower;
        upper[e] = info.upper;
        return e;
    }
    EdgeInfo get_edge_info(Edge e) const {
        return {
            .u = g.u(e), .v=g.v(e),
            .u_head = u_head[e], .v_head = v_head[e],
            .cost_function = cost_function[e],
            .lower = lower[e],
            .upper = upper[e],
        };
    }

public:
    CostScalar cost(Edge e, double x) const;
    CostScalar cost(Solution const&) const;
    /// hack: compute cost of *half* the solution
    CostScalar cost_half(Solution const&) const;
    inline TargetScalar guess(Edge e) const {
        return CostFunction::get_guess(cost_function[e]);
    }
};

struct BiMDFResult {
    std::unique_ptr<BiMDF::Solution> solution;
    BiMDF::CostScalar cost = std::numeric_limits<BiMDF::CostScalar>::quiet_NaN();
};


} // namespace Satsuma
