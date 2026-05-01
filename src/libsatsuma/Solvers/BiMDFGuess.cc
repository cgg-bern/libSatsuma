//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFGuess.hh>
#include <libsatsuma/Problems/CostFunction.hh>

namespace Satsuma {

std::unique_ptr<BiMDF::Guess> make_guess(const BiMDF &bimdf)
{
    auto guessp = std::make_unique<BiMDF::Guess>(bimdf.g);
    auto &guess = *guessp;

    for (const auto e: bimdf.g.edges())
    {
        guess[e] = CostFunction::find_optimal_integer(bimdf.cost_function[e], bimdf.lower[e], bimdf.upper[e]);
    }
    return guessp;
}

} // namespace Satsuma
