//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once
/// support for result serialization using nlohmann-json

#include <nlohmann/json.hpp>
#include <libsatsuma/Extra/Highlevel.hh>
#include <libsatsuma/Solvers/MatchingSolvers.hh>
#include <libTimekeeper/json.hh>

namespace Satsuma {

NLOHMANN_JSON_SERIALIZE_ENUM(BiMCF_to_MCF::Method, {
    {BiMCF_to_MCF::Method::HalfSymmetric, "HalfSymmetric"},
    {BiMCF_to_MCF::Method::HalfAsymmetric, "HalfAsymmetric"},
    {BiMCF_to_MCF::Method::FullSymmetric, "FullSymmetric"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(EveningMode, {
    {EveningMode::MST, "MST"},
    {EveningMode::RoundToEven, "RoundToEven"}
})

NLOHMANN_JSON_SERIALIZE_ENUM(MatchingSolver, {
    {MatchingSolver::Lemon, "Lemon"},
    {MatchingSolver::BlossomV, "BlossomV"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(DeviationLimitKind, {
    {DeviationLimitKind::EdgeFlow, "EdgeFlow"},
    {DeviationLimitKind::NodeThroughflow, "NodeThroughflow"},
})


NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFDoubleCoverInfo,
                                   evening_cost,
                                   evening_n_adjustments,
                                   evening_n_bound_adjustments,
                                   cost,
                                   max_deviation_problem,
                                   max_deviation_solution);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFMatchingInfo,
                                   cost,
                                   cost_changes,
                                   max_refinement_change);


NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFMatchingResult,
                                   // does not store actual solution values!
                                   double_cover_info,
                                   info,
                                   stopwatch);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFperConnectedComponentInfo,
                                   n_nodes,
                                   n_edges,
                                   double_cover,
                                   matching);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFFullResult,
                                   // does not store actual solution values!
                                   cost,
                                   cc_info,
                                   stopwatch);


NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFDoubleCoverConfig,
                                   max_deviation,
                                   matching_solver,
                                   evening_mode,
                                   verbosity,
                                   method);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiMDFSolverConfig,
                                   double_cover,
                                   matching_solver,
                                   refine_with_matching,
                                   refinement_maxdev_min,
                                   refinement_maxdev_max,
                                   deviation_limit,
                                   verbosity);


} // namespace Satsuma
