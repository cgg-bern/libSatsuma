//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "Matching.hh"

namespace Satsuma {


Matching::WeightScalar Matching::cost(Matching::Solution const& _sol) const
{
    Matching::WeightScalar sum = 0;
    for (const auto e: g.edges()) {
        if (_sol[e]) {
            sum += weight[e];
        }
    }
    return sum;
}

bool Matching::is_perfect(const Matching::Solution &sol) const
{
    GraphT::NodeMap<int> count{g, 0};
    for (const auto e: g.edges()) {
        count[g.u(e)] += sol[e];
        count[g.v(e)] += sol[e];
    }
    for (const auto n: g.nodes()) {
        if (count[n] != 1) {
            return false;
        }
    }
    return true;
}

} // namespace Satsuma
