//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMCF.hh>

namespace Satsuma {

BiMCFResult solve_bimcf_gurobi(BiMCF const &bimdf, BiMCF::Guess *x0 = nullptr);

} // namespace Satsuma
