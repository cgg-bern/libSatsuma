//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>

#include <cmath>
#include <cassert>

namespace Satsuma {

BiMDF::CostScalar BiMDF::cost(Edge e, double x) const
{
    return CostFunction::cost(cost_function[e], x);
}

BiMDF::CostScalar BiMDF::cost(Solution const&sol) const
{
    double sum = 0;
    for (auto e: g.edges()) {
        sum += cost(e, sol[e]);
    }
    return sum;
}
BiMDF::CostScalar BiMDF::cost_half(Solution const&sol) const
{
    double sum = 0;
    for (auto e: g.edges()) {
        sum += cost(e, .5 * sol[e]);
    }
    return sum;
}

} // namespace Satsuma
