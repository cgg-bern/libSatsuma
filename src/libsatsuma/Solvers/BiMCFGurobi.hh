#pragma once

#include <libsatsuma/Problems/BiMCF.hh>

namespace Satsuma {

BiMCFResult solve_bimcf_gurobi(BiMCF const &bimdf, BiMCF::Guess *x0 = nullptr);

} // namespace Satsuma
