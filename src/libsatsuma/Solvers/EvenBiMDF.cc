//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Problems/TJoin.hh>
#include <libsatsuma/Solvers/TJoinMST.hh>
#include <libsatsuma/Exceptions.hh>
#include <lemon/maps.h>

namespace Satsuma {

/// return cheapest +/- 1 adjusted target length and the cost of the change
static
std::pair<BiMDF::FlowScalar, BiMDF::CostScalar>
find_best_adjustment(BiMDF const&_bimdf, BiMDF::Edge e, BiMDF::TargetScalar initial)
{
    const auto lower = _bimdf.lower[e];
    const auto upper = _bimdf.upper[e];
    auto guess = std::llround(initial);
    if (guess < lower) {
        guess = lower;
    }
    if (guess > upper) {
        guess = upper;
    }
    const auto base_cost = _bimdf.cost(e, guess);
    auto best_cost = std::numeric_limits<double>::infinity();
    auto best_guess = guess;
    for (int adj = -1; adj<=1; adj+=2)
    {
        int x = guess + adj;

        if (x < lower || x > upper) {
            continue;
        }
        auto cost_change = _bimdf.cost(e, x) - base_cost;
        if (cost_change < best_cost) {
            best_cost = cost_change;
            best_guess = x;
        }
    }
    if (best_guess == guess) {
        throw InternalError("Could not find adjustment for DC; probably upper=lower, not handled yet. lower = " 
                + std::to_string(lower)
                + ", upper = " + std::to_string(upper)
                + ", guess = " + std::to_string(guess));
    }
    assert(best_guess != guess); // bounds too tight
    return {best_guess, best_cost};
}

EveningResult round_to_even(const BiMDF &bimdf)
{
    auto guessp = std::make_unique<BiMDF::Guess>(bimdf.g);

    size_t n_adjustments = 0;
    size_t n_bound_adjustments = 0;
    double cost = 0;
    for (const auto e: bimdf.g.edges())
    {
        const auto lower = bimdf.lower[e];
        const auto upper = bimdf.upper[e];
        BiMDF::FlowScalar guess = std::llround(bimdf.guess(e));
        if (guess < lower) {
            guess = lower;
        } else if (guess > upper) {
            guess = upper;
        }

        if (guess & 1) {
            auto adj = find_best_adjustment(bimdf, e, guess);
            guess = adj.first;
            cost += adj.second;
            ++n_adjustments;
        }
        assert((guess & 1) == 0);
        assert(guess <= bimdf.upper[e]);
        assert(guess >= bimdf.lower[e]);
        (*guessp)[e] = guess;
    }

    return { .guess = std::move(guessp),
             .cost = cost,
             .n_adjustments = n_adjustments,
             .n_bound_adjustments = 0 // TODO?
    };
}



EveningResult TJoinBasedRounding::solve()
{
    const auto &g = bimdf_.g;
    auto guessp = std::make_unique<BiMDF::Guess>(g);
    auto &guess = *guessp;


    for (const auto e: g.edges())
    {
        //assert(bimdf_.upper[e] == BiMDF::inf());
        if (bimdf_.upper[e] <= bimdf_.lower[e]) {
            throw std::runtime_error("BiMDF has arc with upper<= lower, please remove.");
        }
        // TODO: we might round towards an odd bound, which is later rounded towards guess! is this a problem?
        // TODO: use find_best_adjustment
        auto lower = bimdf_.lower[e];
        auto upper = bimdf_.upper[e];
        auto opti = std::llround(bimdf_.guess(e));
        auto rounded = opti;
        if (rounded < lower) {
            rounded = lower;
        }
        guess[e] = rounded;
        auto base_cost = bimdf_.cost(e, rounded);
        cost_[e] = std::numeric_limits<double>::infinity();
        for (int adj = -1; adj<=1; adj+=2)
        {
            int x = rounded + adj;

            if (x < lower || x > upper) {
                continue;
            }
            auto cost_change = bimdf_.cost(e, x) - base_cost;
            if (cost_change < cost_[e]) {
                cost_[e] = cost_change;
                adjusted_guess_[e] = x;
            }
        }
#if 0
        // TODO: try idea: lower cost if guess-lower is odd, we'd like to fix this
        if (std::abs(guess[e] - bimdf_.lower[e]) & 1) {
            cost[e] -= bimdf_.weight[e];
        }
#endif

        // apply flow:
        rhs_[g.u(e)] += bimdf_.u_head[e] ? rounded : -rounded;
        rhs_[g.v(e)] += bimdf_.v_head[e] ? rounded : -rounded;
    }

    for (const auto n: g.nodes())
    {
        tjoin_t_[n] = rhs_[n] & 1;
    }
    auto res = solve_tjoin_mst(tjoin_);

    auto const &sol = *res.solution;

    double cost = 0;
    size_t n_adjustments = 0;
    for (const auto e: g.edges()) {
        if (sol[e]) {
            guess[e] = adjusted_guess_[e];
            auto opti = std::llround(bimdf_.guess(e));
            cost += bimdf_.cost(e, guess[e]) - bimdf_.cost(e, opti);
            ++n_adjustments;
        }
    }


    if (verbosity_ >= 2) {
        std::cerr << "performed " << n_adjustments << " adjustments for even rhs" << std::endl;
    }

#if 0 // rhs is not touched, need to recompute it for this chjeck:
    for (const auto n: g.nodes()) {
        assert((rhs_[n] & 1) == 0);
    }
#endif

#if 0
    // adjust `lower` such that `target-lower` is even
    size_t n_lower_changes = 0;
    for (const auto e: g.edges())
    {
        auto diff = guess[e] - result.lower[e];
        if (diff & 1) {
            ++result.lower[e];
            ++n_lower_changes;
        }
    }

    if (verbosity > 0) {
        std::cerr << "incremented " << n_lower_changes << " lower bounds for even guess-lower" << std::endl;
    }
#endif

    return { .guess = std::move(guessp),
             .cost = cost,
             .n_adjustments = n_adjustments,
             .n_bound_adjustments = 0};
}

EveningResult guess_for_even_rhs(const BiMDF &bimdf, int verbosity, EveningMode mode)
{
    if (mode == EveningMode::RoundToEven) {
      return round_to_even(bimdf);
    } else if (mode == EveningMode::MST) {
        return TJoinBasedRounding(bimdf, verbosity).solve();
    } else {
        assert(false);
        return {};
    }
}



} // namespace Satsuma
