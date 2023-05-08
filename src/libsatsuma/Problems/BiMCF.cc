#include "BiMCF.hh"

namespace Satsuma {

BiMCF::CostScalar BiMCF::compute_cost(const Solution &sol) const
{
    BiMCF::CostScalar sum = 0;
    for (auto e: g.edges()) {
        sum += cost[e] * sol[e];
    }
    return sum;
}


} // namespace Satsuma

