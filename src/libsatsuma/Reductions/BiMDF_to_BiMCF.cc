//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>

#include <cassert>
#include <cmath>

namespace Satsuma {

int add_capacities(int a, int b) {
  if (a == BiMCF::inf() || b == BiMCF::inf()) {
    return BiMCF::inf();
  }
  return a+b;
}

BiMDF_to_BiMCF::BiMDF_to_BiMCF(const BiMDF &_bimdf, Config const& _config)
    : bimdf_(_bimdf)
    , guess_(bimdf_.g) // TODO: _config.guess if copy constructor exists
    , is_forward_(bimcf_.g)
    , mdf_edge_id_(bimcf_.g)
{
    if (_config.out_orig_node) {
        *_config.out_orig_node = std::make_unique<BiMCF::NodeMap<BiMDF::Node>>(bimcf_.g);
    }
    const auto n_nodes = bimdf_.g.maxNodeId() + 1;
    const auto n_edges = bimdf_.g.maxEdgeId() + 1;
#if 0
    // TODO: update the sizes to reserve
    bimcf_.g.reserveNode(n_nodes);
    switch(objective)
    {
    case Objective::SimpleAbs:
        bimcf_.g.reserveEdge(n_edges * 2);
        break;
    case Objective::ExactAbsRefinement:
    case Objective::ExactQuadraticRefinement:
        assert(max_deviation < BiFlowGraph::inf());
        bimcf_.g.reserveEdge(n_edges * 2 * max_deviation);
        break;
    }
#endif


    for (auto &n: bimdf_.g.nodes())
    {
        auto mcf_node = bimcf_.g.addNode();
        if (_config.out_orig_node) {
            (**_config.out_orig_node)[mcf_node] = mcf_node; // not a bug: nodes are identity-mapped, but added in reverse order with ListGraph
        }
        bimcf_.demand[mcf_node] = 0;
    }

    double last_arc_cost = std::numeric_limits<double>::quiet_NaN();
    bool last_forward = false;
    BiMDF::Edge last_mdf_edge = lemon::INVALID;
    BiMDF::Edge last_edge = lemon::INVALID;
    auto add_edge = [&](BiMDF::Edge mdf_edge, bool forward, double cost, int upper = BiMCF::inf())
      -> BiMCF::Edge
    {
        if (upper <= 0)
            return lemon::INVALID;

        if (_config.consolidate &&
            last_mdf_edge == mdf_edge &&
            last_forward == forward &&
            std::fabs(last_arc_cost - cost) <= std::fabs(1e-6 * cost)
            //&& false
            )
        {
          assert(last_edge != lemon::INVALID);
          //std::cout << " -> extending capacity by " << upper;
          bimcf_.upper[last_edge] = add_capacities(bimcf_.upper[last_edge], upper);
          //std::cout << std::endl;
          return last_edge;
        } else {
          auto u = bimdf_.g.u(mdf_edge);
          auto v = bimdf_.g.v(mdf_edge);
          auto e = bimcf_.g.addEdge(u, v);
          mdf_edge_id_[e] = bimdf_.g.id(mdf_edge);
          is_forward_[e] = forward;
          bimcf_.u_head[e] = bimdf_.u_head[mdf_edge] ^ !forward;
          bimcf_.v_head[e] = bimdf_.v_head[mdf_edge] ^ !forward;
          bimcf_.cost[e] = cost;
          bimcf_.upper[e] = upper;

          last_mdf_edge = mdf_edge;
          last_arc_cost = cost;
          last_edge = e;
          last_forward = forward;

          //std::cout << " -> new edge" << std::endl;
          return e;
        }
    };

    for (const auto mdf_edge: bimdf_.g.edges())
    {
        const auto guess = _config.guess[mdf_edge];
        guess_[mdf_edge] = guess;
        auto u = bimdf_.g.u(mdf_edge);
        auto v = bimdf_.g.v(mdf_edge);
        // TODO: applying flow should be a BiFlowGraph member function:
        bimcf_.demand[u] += bimdf_.u_head[mdf_edge] ? -guess : guess;
        bimcf_.demand[v] += bimdf_.v_head[mdf_edge] ? -guess : guess;

        //const double weight = bimdf_.weight[mdf_edge];
        //const double target = bimdf_.target[mdf_edge];

        auto energy = [&](int val) {
            return bimdf_.cost(mdf_edge, val);
        };
        const auto lower = bimdf_.lower[mdf_edge];
        const auto upper = bimdf_.upper[mdf_edge];
        const int cap = _config.even ? 2 : 1;
        const double guess_cost = energy(guess);
        double ecost = guess_cost;

        int dev = 0;


        // forward arcs:
        for (int i = cap; i <= _config.max_deviation; i += cap) {
            int remain = upper - guess - (i - cap); // remaining capacity capacity after applying all *previous* arcs
            int remcap = std::min(cap, remain);
            if (remcap <= 0)
                break;
            auto last_cost = ecost;
            dev = i - cap + remcap;
            ecost = energy(guess+dev);
            double arc_cost = (ecost - last_cost)/remcap;
#if 0
            std::cout << "CCCC cost for forward edge "
                      << bimdf_.g.id(u) << " - " << bimdf_.g.id(v)
                      << ", i = " << i
                      << ", upper = " << upper
                      << ", remcap = " << remcap
                      << ", last_cost = " << last_cost
                      << ", arc_cost = " << arc_cost
                      //<< ", last_arc_cost = " << last_arc_cost
                      << ", ecost =     " << ecost
                      << ", result:" << arc_cost
                      << std::endl;
#endif
            add_edge(mdf_edge, true, arc_cost, remcap);
        }

        const int last_arc_dx = 10;
        if (_config.last_arc_uncapacitated) {
            auto cost = (energy(guess + dev + last_arc_dx)
                       - energy(guess + dev)) / last_arc_dx;
            assert(cost >= 0);
            if (cost < 0 ) { cost = 0;}  // we never want an unbounded arc with negative costs!
            int remain = upper;
            if (upper < BiMCF::inf()) {
                remain -= dev;
            }
            add_edge(mdf_edge, true, cost, remain);
        }

        // backwards arcs:
        ecost = guess_cost;
        for (int i = cap; i <= _config.max_deviation; i += cap) {
            int remain = guess - lower - (i-cap); // remaining capacity capacity after applying all *previous* arcs
            int remcap = std::min(cap, remain);
            if (remcap <= 0)
                break;
            auto last_cost = ecost;
            dev = i - cap + remcap;
            ecost = energy(guess - dev);
            double arc_cost = (ecost - last_cost)/remcap;
#if 0
            std::cout << "DDDD cost for backward edge "
                      << bimdf_.g.id(u) << " - " << bimdf_.g.id(v)
                      << ", i = " << i
                      << ", upper = " << upper
                      << ", remcap = " << remcap
                      << ", last_cost = " << last_cost
                      << ", arc_cost = " << arc_cost
                      //<< ", last_arc_cost = " << last_arc_cost
                      << ", ecost =     " << ecost
                      << ", result:" << arc_cost
                      << std::endl;
#endif
            add_edge(mdf_edge, false, arc_cost, remcap);
        }
        auto remain = guess - lower - dev;
        if (_config.last_arc_uncapacitated && remain > 0) {
            auto cost = (energy(guess - dev - remain)
                       - energy(guess - dev)) / remain;
            add_edge(mdf_edge, false, cost, remain);
        }
    }

#if 0
    std::cerr << "Bi-MDF |V| = " << bimdf_.g.maxNodeId()+1
              << ", |E| = " << bimdf_.g.maxEdgeId()+1
              << std::endl;
    std::cerr << "Bi-MCF |V| = " << bimcf_.g.maxNodeId()+1
              << ", |E| = " << bimcf_.g.maxEdgeId()+1
              << std::endl;
#endif
}

BiMDFResult BiMDF_to_BiMCF::translate_solution(const BiMCFResult &bimcf_res, bool double_guess) const
{
    auto &bimcf_sol = *bimcf_res.solution;

    auto bimdf_solp = std::make_unique<BiMDF::Solution>(bimdf_.g);
    auto &bimdf_sol = *bimdf_solp;

    for (const auto mdf_edge: bimdf_.g.edges())
    {
        bimdf_sol[mdf_edge] = guess_[mdf_edge];
        if (double_guess) {
            bimdf_sol[mdf_edge] *= 2;
        }
    }
    for (auto mcf_edge: bimcf_.g.edges()) {
        auto mdf_edge = bimdf_.g.edgeFromId(mdf_edge_id_[mcf_edge]);
        auto flow = bimcf_sol[mcf_edge];
        bimdf_sol[mdf_edge] += is_forward_[mcf_edge] ? flow : - flow;
    }

    return {.solution = std::move(bimdf_solp)};
}

} // namespace Satsuma
