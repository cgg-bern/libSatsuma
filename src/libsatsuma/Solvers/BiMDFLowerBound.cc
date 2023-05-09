//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFLowerBound.hh>
#include <libsatsuma/Solvers/BiMDFDoubleCover.hh>
#include <libsatsuma/Solvers/BiMDFGuess.hh>
#include <libsatsuma/Solvers/MCF.hh>

namespace Satsuma {
BiMDFLowerBoundResult bimdf_lower_bound(const BiMDF &bimdf,
                                        int initial_maxdev)
{
    auto guessp = make_guess(bimdf);

    assert(initial_maxdev > 1);
    BiMDF::FlowScalar max_dev = initial_maxdev;

    auto best_cost = std::numeric_limits<MCF::CostScalar>::max();
    while (true) {
        auto red_bimcf = BiMDF_to_BiMCF(bimdf, {
                .guess = *guessp,
                .max_deviation = max_dev,
                .last_arc_uncapacitated = true,
                .even = false,
                .consolidate = true});
        auto red_mcf = BiMCF_to_MCF(red_bimcf.bimcf(), {
                                        .method = BiMCF_to_MCF::Method::NotEven});

        MCFResult sol_mcf;

        sol_mcf = solve_mcf_via_lemon_netsimp(red_mcf.mcf());
        auto sol_bimcf = red_mcf.translate_solution(sol_mcf);
        auto sol_bimdf = red_bimcf.translate_solution(sol_bimcf, true);
        // due to integer rounding, the bi-mcf cost may oscillate, so compare integer mcf cost.
        auto cost = sol_mcf.cost;
#if 0
        std::cout << "LB " << sol_bimcf.max_flow
                  << " / " << max_dev
                  << ": " << cost
                  << ", bi-mcf cost " << sol_bimcf.cost
                  << " or " << red_bimcf.bimcf().compute_cost(*sol_bimcf.solution)
                  << ", mcf cost " << sol_mcf.cost
                  << " or " << red_mcf.mcf().compute_cost(*sol_mcf.solution)
                  //<< ", full " << bimdf.cost(*sol_bimdf.solution)
                  << std::endl;
#endif
        if (cost != best_cost) {
            best_cost = cost;
            max_dev *= 2;
        } else {
            return {
                .cost = bimdf.cost_half(*sol_bimdf.solution),
                .max_deviation = sol_bimcf.max_flow
            };
        }
    }
}

} // namespace Satsuma
