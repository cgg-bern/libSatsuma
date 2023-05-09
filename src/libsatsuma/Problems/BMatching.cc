//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "BMatching.hh"
#include "Matching.hh"


namespace Satsuma {

bool BMatching::is_valid(const BMatching::Solution &sol) const
{
    GraphT::NodeMap<int> sum{g, 0};
    for (const auto e: g.edges()) {
        sum[g.u(e)] += sol[e];
        sum[g.v(e)] += sol[e];
    }
    for (const auto n: g.nodes()) {
        if (degree[n] != sum[n]) {
            return false;
        }
    }
    return true;
}

} // namespace Satsuma
