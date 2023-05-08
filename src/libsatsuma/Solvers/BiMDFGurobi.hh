#pragma once

#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {

BiMDFResult solve_bimdf_gurobi(BiMDF const &bimdf, bool int_targets=false, BiMDF::Guess *x0 = nullptr);

BiMDFResult solve_bimdf_gurobi_mcf(BiMDF const &bimdf, BiMDF::Guess *x0 = nullptr);

} // namespace Satsuma
