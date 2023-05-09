//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Solvers/BiMDFDoubleCover.hh>
#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Solvers/Matching.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#if SATSUMA_HAVE_GUROBI
#include <libsatsuma/Config/Gurobi.hh>
#endif
#include <libTimekeeper/StopWatch.hh>

namespace Satsuma {

struct BiMDFMatchingInfo {
    BiMDF::CostScalar cost;
    std::vector<double> cost_changes;
    int max_refinement_change;
};

struct BiMDFMatchingResult {
    BiMDFResult result;
    BiMDFDoubleCoverInfo double_cover_info;
    BiMDFMatchingInfo info;
    Timekeeper::HierarchicalStopWatchResult stopwatch;
};

struct BiMDFSolverConfig {
    BiMDFDoubleCoverConfig double_cover;
    /// matching solver to use for refinement (can theoretically be different from solver used for DC)
    MatchingSolver matching_solver = MatchingSolver::Default;

    bool refine_with_matching = true;
    /// Maximum deviation from x0 in primary iterations.
    /// 1 or 2 are recommended.
    int refinement_maxdev_min = 2;
    /// Maximum deviation from x0 in last iteration.
    /// Note: 2 always suffices for an exact solution.
    int refinement_maxdev_max = 2;
    DeviationLimitKind deviation_limit = DeviationLimitKind::Default;
    int verbosity = 2;
};


BiMDFMatchingResult solve_bimdf_matching(
        const BiMDF &bimdf,
        BiMDFSolverConfig const& _config = BiMDFSolverConfig());

struct BiMDFperConnectedComponentInfo {
    size_t n_nodes;
    size_t n_edges;
    BiMDFDoubleCoverInfo double_cover;
    BiMDFMatchingInfo matching;
};

struct BiMDFFullResult {
    std::unique_ptr<BiMDF::Solution> solution;
    BiMDF::CostScalar cost;
    std::vector<BiMDFperConnectedComponentInfo> cc_info;
    Timekeeper::HierarchicalStopWatchResult stopwatch;

};

/// Solve BiMDF by solving each connected component separately, applying
/// simplification on each one
BiMDFFullResult solve_bimdf(const Satsuma::BiMDF &_bimdf,
                            BiMDFSolverConfig const &_config = BiMDFSolverConfig());


} // namespace Satsuma
