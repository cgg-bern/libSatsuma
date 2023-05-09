//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Solvers/Matching.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#if SATSUMA_HAVE_GUROBI
#include <libsatsuma/Config/Gurobi.hh>
#endif
#include <libTimekeeper/StopWatch.hh>

namespace Satsuma {

struct BiMDFDoubleCoverInfo {
    TJoin::CostScalar evening_cost;
    size_t evening_n_adjustments;
    size_t evening_n_bound_adjustments;
    BiMDF::CostScalar cost;
    BiMCF::FlowScalar max_deviation_problem;
    BiMCF::FlowScalar max_deviation_solution;
};

struct BiMDFDoubleCoverResult {
    std::unique_ptr<BiMDF::Solution> solution;
    BiMDFDoubleCoverInfo info;
    Timekeeper::HierarchicalStopWatchResult stopwatch;
};

struct BiMDFDoubleCoverConfig {
    /// Maximum deviation from initial guess that is represented exactly
    int max_deviation = 5;
    /// Used in MST evening_mode
    MatchingSolver matching_solver = MatchingSolver::Default;
    EveningMode evening_mode = EveningMode::MST;
    int verbosity = 2;
    BiMCF_to_MCF::Method method = BiMCF_to_MCF::Method::Default;
};


BiMDFDoubleCoverResult approximate_bimdf_doublecover(
        const BiMDF &bimdf,
        BiMDFDoubleCoverConfig const& config = BiMDFDoubleCoverConfig());

} // namespace Satsuma
