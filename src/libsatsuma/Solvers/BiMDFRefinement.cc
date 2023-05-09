//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFRefinement.hh>
#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Solvers/Matching.hh>
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

BiMDFRefinementResult refine_with_matching(const BiMDF &_bimdf,
                                           BiMDF::Solution const& f0,
                                           int max_deviation,
                                           DeviationLimitKind deviation_limit,
                                           MatchingSolver matching_solver)
{
    auto red_bimcf = BiMDF_to_BiMCF(_bimdf, {
        .guess = f0,
        .max_deviation = max_deviation,
        .last_arc_uncapacitated = false, // matching is not compatible with uncapacitated arcs
        .even = false,
        .consolidate = true});
#if 0
    for (auto n: red_bimcf.bimcf().g.nodes()) {
        if(red_bimcf.bimcf().demand[n] != 0) {
            throw std::runtime_error("f0 must be conservative");
        }
    }
#endif

#if 0
    auto sol_bimcf = solve_bimcf_gurobi(red_bimcf.bimcf());
#else
    auto red_bmatching = BiMCF_to_BMatching(red_bimcf.bimcf(), {
            .max_deviation = max_deviation,
            .deviation_limit = deviation_limit});
    auto red_matching = BMatching_to_Matching(red_bmatching.bmatching());

    auto sol_matching = solve_matching(red_matching.matching(), matching_solver);
    auto sol_bmatching = red_matching.translate_solution(sol_matching);
    auto sol_bimcf = red_bmatching.translate_solution(sol_bmatching);
#endif
    auto sol_bimdf = red_bimcf.translate_solution(sol_bimcf);
    return {.sol = std::move(sol_bimdf.solution),
            .cost_change = sol_bimcf.cost};

}


} // namespace Satsuma
