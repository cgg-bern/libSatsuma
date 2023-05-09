//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFDoubleCover.hh>
#include <libsatsuma/Solvers/BiMDFGuess.hh>
#include <libsatsuma/Solvers/MCF.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#include <libsatsuma/Reductions/BMatching_to_Matching.hh>
#include <libsatsuma/Reductions/BiMDF_Simplification.hh>
#include <libsatsuma/Exceptions.hh>

#if SATSUMA_HAVE_GUROBI
#  include <libsatsuma/Solvers/BiMCFGurobi.hh> // just for testing
#endif

namespace Satsuma {


BiMDFDoubleCoverResult
approximate_bimdf_doublecover(
        const BiMDF &_bimdf,
        const BiMDFDoubleCoverConfig &_config)
{
    Timekeeper::HierarchicalStopWatch sw_root{"bimdf via DC"};
    Timekeeper::HierarchicalStopWatch sw_evening{"evening", sw_root};
    Timekeeper::HierarchicalStopWatch sw_reductions{"reductions", sw_root};
    Timekeeper::HierarchicalStopWatch sw_solve{"solve", sw_root};

    sw_root.resume();
    sw_evening.resume();
    auto evening = guess_for_even_rhs(_bimdf, _config.verbosity, _config.evening_mode);
    sw_evening.stop();

    sw_reductions.resume();
    auto red_bimcf = BiMDF_to_BiMCF(_bimdf, {
            .guess = *evening.guess,
            .max_deviation = _config.max_deviation,
            .last_arc_uncapacitated = true,
            .even = true,
            .consolidate = true});
    auto red_mcf = BiMCF_to_MCF(red_bimcf.bimcf(), {
                                    .method = _config.method});
    sw_reductions.stop();

    sw_solve.resume();
    auto sol_mcf = solve_mcf_via_lemon_netsimp(red_mcf.mcf());
    sw_solve.stop();

    sw_reductions.resume();
    auto sol_bimcf = red_mcf.translate_solution(sol_mcf);
#if 0 // not true when we have zero-cost cycles
    if (verbosity > 1 && sol_bimcf.max_flow >= max_deviation) {
        std::cerr << "WARNING: DC appromixation: max_deviation too low, cost function not represented exactly: "
                  << sol_bimcf.max_flow << " / " << max_deviation
                  << std::endl;
    }
#endif
    //assert(red_bimcf.bimcf().is_valid(*sol_bimcf.solution));
    auto sol_bimdf = red_bimcf.translate_solution(sol_bimcf);
    //assert(bimdf.is_valid(*sol_bimdf.solution));
    sw_reductions.stop();
    if (!_bimdf.is_valid(*sol_bimdf.solution)) {
        throw InternalError("approximation result infeasible.");
    }

    auto cost = _bimdf.cost(*sol_bimdf.solution);
    sw_root.stop();
    return {.solution = std::move(sol_bimdf.solution),
            .info = {
                .evening_cost = evening.cost,
                .evening_n_adjustments = evening.n_adjustments,
                .evening_n_bound_adjustments = evening.n_bound_adjustments,
                .cost = cost,
                .max_deviation_solution = sol_bimcf.max_flow,
            },
            .stopwatch = sw_root};
}

} // namespace Satsuma
