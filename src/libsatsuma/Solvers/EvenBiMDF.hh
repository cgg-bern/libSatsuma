//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/TJoin.hh>

#include <cmath>
#include <cassert>
#include <memory>


namespace Satsuma {

enum class EveningMode {
  MST,        // T-Join based rounding with MST approximation
  RoundToEven // round each target length to an even value. Requires that odd lengths can always be adjusted (lower<upper)
};

struct EveningResult {
    std::unique_ptr<BiMDF::Guess> guess;
    double cost;
    size_t n_adjustments;
    size_t n_bound_adjustments;
};


EveningResult round_to_even(const BiMDF &mdf);

class TJoinBasedRounding {
public:
    TJoinBasedRounding(const BiMDF &_bimdf, int _verbosity=2)
        : bimdf_(_bimdf)
        , verbosity_(_verbosity)
        , cost_(bimdf_.g)
        , rhs_(bimdf_.g, 0)
        , adjusted_guess_(bimdf_.g)
        , tjoin_t_(bimdf_.g)
        , tjoin_(bimdf_.g, tjoin_t_, cost_)
    {}

    EveningResult solve();

    TJoin const& tjoin() const {return tjoin_;}

private:
    BiMDF const& bimdf_;
    int verbosity_;
    BiMDF::EdgeMap<TJoin::CostScalar> cost_;
    BiMDF::NodeMap<int> rhs_;
    BiMDF::EdgeMap<int> adjusted_guess_;
    TJoin::TMap tjoin_t_;
    TJoin tjoin_;

};
EveningResult guess_for_even_rhs(const BiMDF &bimdf, int verbosity, EveningMode mode);


} // namespace Satsuma

