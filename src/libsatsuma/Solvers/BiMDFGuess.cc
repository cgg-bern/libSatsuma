//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFGuess.hh>

namespace Satsuma {

std::unique_ptr<BiMDF::Guess> make_guess(const BiMDF &bimdf)
{
    auto guessp = std::make_unique<BiMDF::Guess>(bimdf.g);
    auto &guess = *guessp;

    for (const auto e: bimdf.g.edges())
    {
        auto target = bimdf.guess(e);

        auto floor = static_cast<BiMDF::FlowScalar>(std::llround(std::floor(target)));
        auto ceil = static_cast<BiMDF::FlowScalar>(std::llround(std::ceil(target)));
        floor = std::max(floor, bimdf.lower[e]);
        ceil = std::min(floor, bimdf.upper[e]);
        if(floor > ceil) {
            guess[e] = floor;
        } else {
            auto min_cost = std::numeric_limits<BiMDF::CostScalar>::infinity();
            auto last_cost = std::numeric_limits<BiMDF::CostScalar>::infinity();
            for (auto x = floor; x <= ceil; ++x) {
                auto cost = bimdf.cost(e, x);
                if (cost > last_cost) {
                    // increasing cost -> we are past the minimum of the convex cost function
                    break;
                }
                last_cost = cost;
                if (cost < min_cost) {
                    min_cost = cost;
                    guess[e] = x;
                }
            }
        }
    }
    return guessp;
}

} // namespace Satsuma
