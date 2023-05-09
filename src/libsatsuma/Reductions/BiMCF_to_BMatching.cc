//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#include <cassert>
#include <stdexcept>

namespace Satsuma {

BiMCF_to_BMatching::BiMCF_to_BMatching(const BiMCF &_bimcf,
                                       const Config &_config)
    : bimcf_(_bimcf)
{

    if (_config.out_orig_node) {
        *_config.out_orig_node = std::make_unique<BMatching::NodeMap<BiMCF::Node>>(bmatching_.g);
    }
    if (_config.out_is_in_node) {
        *_config.out_is_in_node = std::make_unique<BMatching::NodeMap<bool>>(bmatching_.g);
    }
    if (_config.out_orig_edge) {
        *_config.out_orig_edge = std::make_unique<BMatching::EdgeMap<BiMCF::Edge>>(bmatching_.g);
    }


    auto max_deviation = _config.max_deviation;
    auto node_in = [&](BiMCF::Node const&n) -> BMatching::Node {
        return bmatching_.g.nodeFromId(bimcf_.g.id(n) * 2);
    };
    auto node_out = [&](BiMCF::Node const&n) -> BMatching::Node {
        return bmatching_.g.nodeFromId(bimcf_.g.id(n) * 2 + 1);
    };

    size_t n_bounded_edges = 0;
    size_t n_unbounded_edges = 0;
    for (auto e: bimcf_.g.edges()) {
        if (bimcf_.upper[e] < BiMCF::inf()) {
            ++n_bounded_edges;
        } else {
            ++n_unbounded_edges;
        }
    }

    const auto n_nodes = _bimcf.g.maxNodeId() + 1;
    bmatching_.g.reserveNode(n_nodes * 2 + n_bounded_edges * 2);
    bmatching_.g.reserveEdge(n_unbounded_edges + 3 * n_bounded_edges);

    for ([[maybe_unused]] const auto bimcf_node: _bimcf.g.nodes())
    {
        bmatching_.g.addNode();
        bmatching_.g.addNode();
    }
    // how much additional flow can enter or leave each node?
    BMatching::NodeMap<int> max_flow_in(_bimcf.g, 0);
    BMatching::NodeMap<int> max_flow_out(_bimcf.g, 0);
    if (_config.deviation_limit == DeviationLimitKind::EdgeFlow)
    {
        for (const auto e: _bimcf.g.edges())
        {
            auto max_inc = std::min(_bimcf.upper[e], _config.max_deviation);
            if (_bimcf.u_head[e]) {
                max_flow_in[_bimcf.g.u(e)] += max_inc;
            } else {
                max_flow_out[_bimcf.g.u(e)] += max_inc;
            }
            if (_bimcf.v_head[e]) {
                max_flow_in[_bimcf.g.v(e)] += max_inc;
            } else {
                max_flow_out[_bimcf.g.v(e)] += max_inc;
            }
        }
    }
    for (const auto bimcf_node: _bimcf.g.nodes())
    {
        auto n_in = node_in(bimcf_node);
        auto n_out = node_out(bimcf_node);
        auto max_node_flow = max_deviation;
        if (_config.deviation_limit == DeviationLimitKind::EdgeFlow)
        {
            // this assumes zero demand
            max_node_flow = std::min(
                    max_flow_in[bimcf_node],
                    max_flow_out[bimcf_node]);
        }

        auto demand = bimcf_.demand[bimcf_node];
        //assert(demand == 0); // non-zero demands not tested and likely broken
        if (demand != 0) {
            throw std::runtime_error("BiMCF_to_BMatching: non-zero demands currently not supported");
        }
        if (demand > 0) {
            bmatching_.degree[n_in] = max_node_flow + demand;
            bmatching_.degree[n_out] = max_node_flow;
        } else {
            bmatching_.degree[n_in] = max_node_flow;
            bmatching_.degree[n_out] = max_node_flow - demand;
        }
        auto e = bmatching_.g.addEdge(n_in, n_out);
        bmatching_.capacity[e] = max_node_flow;
        bmatching_.weight[e] = 0;


        if (_config.out_is_in_node) {
            (**_config.out_is_in_node)[n_in] = true;
            (**_config.out_is_in_node)[n_out] = false;
        }
        if (_config.out_orig_node) {
            (**_config.out_orig_node)[n_in] = bimcf_node;
            (**_config.out_orig_node)[n_out] = bimcf_node;
        }
        if (_config.out_orig_edge) {
            (**_config.out_orig_edge)[e] = lemon::INVALID;
        }
    }

    for (auto mcf_edge: bimcf_.g.edges()) {
        auto u = bimcf_.u_head[mcf_edge] ? node_in(bimcf_.g.u(mcf_edge)) : node_out(bimcf_.g.u(mcf_edge));
        auto v = bimcf_.v_head[mcf_edge] ? node_in(bimcf_.g.v(mcf_edge)) : node_out(bimcf_.g.v(mcf_edge));
        auto upper = bimcf_.upper[mcf_edge];
        {
            // not bounded, add regular edge
            auto e = bmatching_.g.addEdge(u, v);
            bmatching_.weight[e] = -bimcf_.cost[mcf_edge];
            bmatching_.capacity[e] = upper;
            bm_edge_id[mcf_edge] = bmatching_.g.id(e);

            if (_config.out_orig_edge) {
                (**_config.out_orig_edge)[e] = mcf_edge;
            }
        }
    }
}

BiMCFResult BiMCF_to_BMatching::translate_solution(const BMatchingResult &bmatching_result) const
{
    const auto &bsol = *bmatching_result.solution;

    auto sol = std::make_unique<BiMCF::Solution>(bimcf_.g, 0);
    for (auto mcf_edge: bimcf_.g.edges()) {
        auto bm_edge = bmatching_.g.edgeFromId(bm_edge_id[mcf_edge]);
        (*sol)[mcf_edge] = bsol[bm_edge];
    }
    if (!bimcf_.is_valid(*sol)) {
        throw std::runtime_error("mcf solution kaput");
    }
    return {.solution = std::move(sol), .cost = -bmatching_result.weight};


}

} // namespace Satsuma
