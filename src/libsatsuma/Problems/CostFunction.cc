#include <libsatsuma/Problems/CostFunction.hh>
#include <algorithm>
#include <cassert>

namespace Satsuma::CostFunction {

BaseObjective::~BaseObjective() = default;

int find_optimal_integer(CostFunction::Function const& _f, int _lower, int _upper)
{
    using FlowScalar = int;
    using CostScalar = double;
    auto target = CostFunction::get_guess(_f);
    auto floor = static_cast<FlowScalar>(std::llround(std::floor(target)));
    auto ceil = static_cast<FlowScalar>(std::llround(std::ceil(target)));
    floor = std::max(floor, _lower);
    ceil = std::min(floor, _upper);
    FlowScalar guess = floor;
    if(floor > ceil) {
        assert(false);
        // TODO: emit warning!
        return floor;
    } else {
        auto min_cost = std::numeric_limits<CostScalar>::infinity();
        auto last_cost = std::numeric_limits<CostScalar>::infinity();
        for (auto x = floor; x <= ceil; ++x) {
            auto cost = CostFunction::cost(_f, x);
            if (cost > last_cost) {
                // increasing cost -> we are past the minimum of the convex cost function
                break;
            }
            last_cost = cost;
            if (cost < min_cost) {
                min_cost = cost;
                guess = x;
            }
        }
    }
    return guess;
}
} // namespace
