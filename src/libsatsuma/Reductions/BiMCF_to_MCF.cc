//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <cmath>
#include <cassert>

namespace Satsuma {

BiMCF_to_MCF::BiMCF_to_MCF(const BiMCF &_bimcf,
                           Config const &_config)
    : bimcf_(_bimcf)
    , method_(_config.method)
{
    static_assert(std::numeric_limits<MCF::CostScalar>::max() >= (1LL<<63)); // for 32-bit ints, adjust costmul
    static_assert(std::is_signed_v<MCF::CostScalar>);
    costmul_ = 1LL << 20; // Warning: if we choose this dynamically, bimdf_lower_bound does not work anymore!
    BiMCF::CostScalar abs_cost_sum = 0;
    for (const auto e: _bimcf.g.edges()) {
        abs_cost_sum += std::fabs(_bimcf.cost[e]);
    }
    if (abs_cost_sum > 1LL<<40) { // conservative guess
        throw std::runtime_error("BiMCF_to_MCF: costs too high");
    }

    if (_config.out_orig_node) {
        *_config.out_orig_node = std::make_unique<MCF::NodeMap<BiMCF::Node>>(mcf_.g);
    }
    if (_config.out_node_is_plus) {
        *_config.out_node_is_plus = std::make_unique<MCF::NodeMap<bool>>(mcf_.g);
    }
    auto node_plus = [&](BiMCF::Node const&n) -> MCF::Node {
        return mcf_.g.nodeFromId(bimcf_.g.id(n) * 2);
    };
    auto node_minus = [&](BiMCF::Node const&n) -> MCF::Node {
        return mcf_.g.nodeFromId(bimcf_.g.id(n) * 2 + 1);
    };

    mcf_.g.reserveNode((bimcf_.g.maxNodeId()+1)*2);
    mcf_.g.reserveArc((bimcf_.g.maxEdgeId()+1)*2);

    for ([[maybe_unused]] auto bimcf_node: bimcf_.g.nodes())
    {
        mcf_.g.addNode();
        mcf_.g.addNode();
    }
    for (auto bimcf_node: bimcf_.g.nodes())
    {
        auto n_plus = node_plus(bimcf_node);
        auto n_minus = node_minus(bimcf_node);

        auto demand = bimcf_.demand[bimcf_node];

        if (method_ != Method::NotEven) {
            assert((demand&1) == 0);
        }
        if (method_ == Method::HalfSymmetric
                  || method_ == Method::HalfAsymmetric)
        {
            demand /= 2;
        }
        mcf_.supply[n_plus] = - demand;
        mcf_.supply[n_minus] = demand;

        if (_config.out_node_is_plus) {
            (**_config.out_node_is_plus)[n_plus] = true;
            (**_config.out_node_is_plus)[n_minus] = false;
        }
        if (_config.out_orig_node) {
            (**_config.out_orig_node)[n_plus] = bimcf_node;
            (**_config.out_orig_node)[n_minus] = bimcf_node;
        }
    }

    for (auto bimcf_edge: bimcf_.g.edges())
    {
        auto mcf_u = bimcf_.g.u(bimcf_edge);
        auto mcf_v = bimcf_.g.v(bimcf_edge);
        bool u_head = bimcf_.u_head[bimcf_edge];
        bool v_head = bimcf_.v_head[bimcf_edge];

        auto u_plus = node_plus(mcf_u);
        auto u_minus = node_minus(mcf_u);
        auto v_plus = node_plus(mcf_v);
        auto v_minus = node_minus(mcf_v);

        assert(bimcf_.lower[bimcf_edge] == 0);

        MCF::Node src0 = u_head ? u_minus : u_plus;
        MCF::Node dst0 = v_head ? v_plus : v_minus;

        MCF::Node src1 = v_head ? v_minus : v_plus;
        MCF::Node dst1 = u_head ? u_plus : u_minus;

#if 0
        if (u_head) {
            if (v_head) {
                // u<->v
                arc0 = {u_minus, v_plus};
                arc1 = {v_minus, u_plus};
            } else {
                // u<--v
                arc0 = {u_minus, v_minus};
                arc1 = {v_plus, u_plus};
            }
        } else {
            if (v_head) {
                // u-->v
                arc0 = {u_plus, v_plus};
                arc1 = {v_minus, u_minus};
            } else {
                // u>-<v
                arc0 = {u_plus, v_minus};
                arc1 = {v_plus, u_minus};
            }
        }
#endif
        double scaled_cost = bimcf_.cost[bimcf_edge] * costmul_;
        auto cost = static_cast<MCF::CostScalar>(std::llround(scaled_cost));
        if (mcf_u == mcf_v) {
          // special case self-loop to avoid double arcs - not necessary, just neat.
          // could be replaced by arc-combining postprocessing step.
          auto a = mcf_.g.addArc(src0, dst0);
          mcf_.cost[a] = cost;
          mcf_.upper[a] = bimcf_.upper[bimcf_edge];
          orig_bimcf_edge_[a] = bimcf_edge;
        } else {

          auto upper = bimcf_.upper[bimcf_edge];
          if (method_ == Method::FullSymmetric
                  || method_ == Method::HalfSymmetric)
          {
              // we strictly require even right hand sides
              if (upper & 1) {
                  --upper;
              }
          }
          auto up0 = upper;
          auto up1 = upper;
          if (upper == BiMCF::inf()) {
              up0 = up1 = BiMCF::inf();
          } else if (method_ == Method::HalfSymmetric
                  || method_ == Method::HalfAsymmetric)
          {
              up0 = upper/2;
              up1 = upper - up0;
          }

          if (up0) {
              auto a = mcf_.g.addArc(src0, dst0);
              mcf_.cost[a] = cost;
              mcf_.upper[a] = up0;
              orig_bimcf_edge_[a] = bimcf_edge;
          }
          if (up1) {
              auto a = mcf_.g.addArc(src1, dst1);
              mcf_.cost[a] = cost;
              mcf_.upper[a] = up1;
              orig_bimcf_edge_[a] = bimcf_edge;
          }
        }
    }

}

BiMCFResult BiMCF_to_MCF::translate_solution(const MCFResult &mcf_result) const
{
    auto sol = std::make_unique<BiMCF::Solution>(bimcf_.g, 0);

    const auto &mcf_sol = *mcf_result.solution;

    // translate solution back:
    BiMCF::CostScalar cost = 0;
    for (auto arc: mcf_.g.arcs()) {
        auto e = orig_bimcf_edge_[arc];
        auto val = mcf_sol[arc];
        if (method_ == Method::FullSymmetric) {
          assert((val & 1) == 0); // TODO: could already halve here, I suppose
        }
        (*sol)[e] += val;
        cost += bimcf_.cost[e] * val;
    }
    BiMCF::FlowScalar max_flow = 0;
    for (auto e: bimcf_.g.edges()) {
        auto &val = (*sol)[e];
        if (method_ == Method::FullSymmetric) {
            assert((val & 1) == 0);
            if (val & 1) {
                throw std::runtime_error("Sum is odd, can't halve! I guess MCF solution is not basic.");
            }
            val /= 2;
        }
        if (val > max_flow) {
            max_flow = val;
        }
    }
    if (method_ == Method::FullSymmetric
            || method_ == Method::NotEven) {
        cost *= .5;
    }
    return {.solution = std::move(sol),
            .cost = cost,
            .max_flow = max_flow};

}

} // namespace Satsuma
