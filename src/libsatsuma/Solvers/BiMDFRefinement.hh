//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Solvers/BiMDFDoubleCover.hh>
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Solvers/Matching.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#include <libTimekeeper/StopWatch.hh>

namespace Satsuma {

struct BiMDFRefinementResult {
    std::unique_ptr<BiMDF::Solution> sol;
    BiMDF::CostScalar cost_change;
};

BiMDFRefinementResult refine_with_matching(const BiMDF &_orig,
                                           BiMDF::Solution const& f0,
                                           int max_change,
                                           DeviationLimitKind deviation_limit = DeviationLimitKind::Default,
                                           MatchingSolver matching_solver = MatchingSolver::Default);

} // namespace Satsuma
