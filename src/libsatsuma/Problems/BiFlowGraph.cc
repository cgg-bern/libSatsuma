//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "BiFlowGraph.hh"
#include <cassert>

namespace Satsuma {

bool BiFlowGraph::is_valid(const BiFlowGraph::Solution &sol) const
{
    GraphT::NodeMap<int> sum{g, 0};
    for (const auto e: g.edges()) {
        if (sol[e] < lower[e]) {
#if DEBUG_INVALID
            std::cout << "lower bound violated on edge "
                      << g.id(e) << ": "
                      << sol[e] << " < " << lower[e] << std::endl;
#endif
            return false;
        }
        if (sol[e] > upper[e]) {
#if DEBUG_INVALID
            std::cout << "upper bound violated on edge "
                      << g.id(e) << ": "
                      << sol[e] << " > " << upper[e] << std::endl;
#endif
            return false;
        }
        sum[g.u(e)] += u_head[e] ? sol[e] : -sol[e];
        sum[g.v(e)] += v_head[e] ? sol[e] : -sol[e];
    }
    for (const auto n: g.nodes()) {
        if (demand[n] != sum[n]) {
#if DEBUG_INVALID
            std::cout << "demand not fulfilled" << std::endl;
#endif
            return false;
        }
    }
    return true;
}

void BiFlowGraph::apply_flow(
        const BiFlowGraph::Guess &f,
        NodeMap<FlowScalar> &out) const
{
    for (const auto n: g.nodes()) {
        out[n] = -demand[n];
    }
    for (const auto e: g.edges()) {
        assert (f[e] >= lower[e]);
        assert (f[e] <= upper[e]);
        out[g.u(e)] += u_head[e] ? f[e] : -f[e];
        out[g.v(e)] += v_head[e] ? f[e] : -f[e];
    }
}
auto BiFlowGraph::apply_flow(const BiFlowGraph::Guess &f) const
    -> std::unique_ptr<NodeMap<FlowScalar>>
{
    auto net = std::make_unique<NodeMap<FlowScalar>>(g);
    apply_flow(f, *net);
    return net;
}


} // namespace Satsuma
