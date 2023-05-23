//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "MCF.hh"


namespace Satsuma {

MCF::CostScalar MCF::compute_cost(const Solution &sol) const
{
    MCF::CostScalar sum = 0;
    for (auto a: g.arcs()) {
        sum += cost[a] * sol[a];
    }
    return sum;
}

MCF::FlowScalar MCF::inflow(MCF::Node n, MCF::Solution const&sol) const
{
    FlowScalar sum = 0;
    for (const auto a: g.inArcs(n)) {
        sum += sol[a];
    }
    return sum;
}

MCF::FlowScalar MCF::outflow(MCF::Node n, MCF::Solution const&sol) const
{
    FlowScalar sum = 0;
    for (const auto a: g.outArcs(n)) {
        sum += sol[a];
    }
    return sum;
}

} // namespace Satsuma
