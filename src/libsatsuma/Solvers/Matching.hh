//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once
#include <libsatsuma/Problems/Matching.hh>
#include <libsatsuma/Solvers/MatchingSolvers.hh>

namespace Satsuma {

/// Use Lemon by default, it is always available
MatchingResult solve_matching(Matching const &mp,
                              MatchingSolver solver = MatchingSolver::Default);

MatchingResult solve_matching_via_lemon(Matching const &mp);

/// The matching problem must be feasible, otherwise Blossom-V warns of undefined behaviour!
MatchingResult solve_matching_via_blossomV(Matching const &mp);

} // namespace Satsuma
