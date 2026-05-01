//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Config/Export.hh>
#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {

SATSUMA_EXPORT BiMDFResult solve_bimdf_gurobi(BiMDF const &bimdf, bool int_targets=false, BiMDF::Guess *x0 = nullptr);

SATSUMA_EXPORT BiMDFResult solve_bimdf_gurobi_mcf(BiMDF const &bimdf, BiMDF::Guess *x0 = nullptr);

} // namespace Satsuma
