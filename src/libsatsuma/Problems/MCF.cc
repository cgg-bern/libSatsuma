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

} // namespace Satsuma
