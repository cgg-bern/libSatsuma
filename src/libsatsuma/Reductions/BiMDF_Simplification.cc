//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Reductions/BiMDF_Simplification.hh>
#include <lemon/connectivity.h>
#include <cassert>

namespace Satsuma {

BiMDF_Simplification::BiMDF_Simplification(const BiMDF &_orig)
    : orig_{_orig}
    , simp_edge_{_orig.g}
{

    BiMDF::NodeMap<bool> collapse_node{_orig.g, false};
    BiMDF::EdgeMap<bool> collapse_edge{_orig.g, false};

    size_t n_collapse_nodes = 0;
    for (const auto n: _orig.g.nodes()) {
        if (_orig.demand[n] != 0) {
            continue;
        }
        bool have_one_head = false;
        bool have_one_tail = false;

        for([[maybe_unused]] const auto e: _orig.g.outArcs(n)) {
            bool is_head = (n == _orig.g.u(e)) ? _orig.u_head[e] : _orig.v_head[e];
            if (is_head) {
                if (have_one_head) {
                    // more than one head
                    have_one_head = false;
                    break;
                }
                have_one_head = true;
            } else {
                if (have_one_tail) {
                    // more than one tail
                    have_one_tail = false;
                    break;
                }
                have_one_tail = true;
            }
        }
        if (have_one_head && have_one_tail) {
            ++n_collapse_nodes;
            collapse_node[n] = true;
            for(const auto e: _orig.g.outArcs(n)) {
                collapse_edge[e] = true;
            }
        }
    }
    std::cout << "# collapsible nodes:" << n_collapse_nodes
              << " / " << _orig.g.maxNodeId() + 1
              << std::endl;

    if (n_collapse_nodes == static_cast<size_t>(_orig.g.maxEdgeId() + 1)) {
        std::cout << "graph is collapsible cycle, preventing collapse for an arbitrary node." << std::endl;
        collapse_node[_orig.g.nodeFromId(_orig.g.maxNodeId())] = false;
    }

    auto simp_node = NodeMap<Node>{_orig.g, lemon::INVALID};
    for (const auto n: _orig.g.nodes()) {
        if (!collapse_node[n]) {

            simp_node[n] = simp_.g.addNode();
        }
    }

    BiMDF::EdgeMap<bool> edge_added{_orig.g, false};
    std::vector<Edge> chain;
    std::vector<CostFunction::Function> chain_costs;
    auto extend_chain = [&](Edge _e, Node _n) -> std::pair<Node, bool> {
        const auto &g = _orig.g;
        Edge cur_e = _e;
        Node cur_n = _n;
        //std::cout << "extend_chain(e = " << g.id(_e) << ", n = " << g.id(_n) << ")" << std::endl;
        while (collapse_node[cur_n]) {
            //std::cout << "\tcur_e = " << g.id(cur_e) << ", cur_n = " << g.id(cur_n) << std::endl;
            for(const auto &other_a: g.outArcs(cur_n)) {
                //std::cout << "\t\tother_e = " << g.id(other_a) << std::endl;
                if (Edge(other_a) == cur_e) {
                    continue;
                } else {
                    chain.push_back(other_a);
                    const auto other_u = g.u(other_a);
                    const auto other_v = g.v(other_a);
                    //std::cout << "\t\t\tother_u = " << g.id(other_u) << std::endl;
                    //std::cout << "\t\t\tother_v = " << g.id(other_v) << std::endl;
                    assert(other_u == cur_n || other_v == cur_n);
                    cur_n = (cur_n == other_u) ? other_v : other_u;

                    cur_e = other_a;
                    if(edge_added[cur_e]) {
                        throw std::runtime_error("this should not happen, found edge twice");
                    }
                    edge_added[cur_e] = true;
                    break; // disable for debug
                }
            }
        }
        bool head = g.u(cur_e) == cur_n ? _orig.u_head[cur_e] : _orig.v_head[cur_e];
        return {cur_n, head};
    };

    for (const auto e: _orig.g.edges()) {
        if (edge_added[e]) continue;
        if (!collapse_edge[e]) {
            BiMDF::EdgeInfo ei = _orig.get_edge_info(e);
            ei.u = simp_node[ei.u];
            ei.v = simp_node[ei.v];
            simp_edge_[e] = simp_.add_edge(ei);
        } else {
            chain.clear();
            chain_costs.clear();
            chain.push_back(e);
            auto [u, u_head] = extend_chain(e, _orig.g.u(e));
            auto [v, v_head] = extend_chain(e, _orig.g.v(e));
            auto ei = BiMDF::EdgeInfo{
                    .u = simp_node[u], .v = simp_node[v],
                    .u_head = u_head, .v_head = v_head};
            if (ei.u == lemon::INVALID || ei.v == lemon::INVALID) {
                throw std::runtime_error("something wrong with extend_chain, end node not in new problem");
            }
            ei.lower = _orig.lower[e];
            ei.upper = _orig.upper[e];
            for (const auto chain_e: chain) {
                ei.lower = std::max(ei.lower, _orig.lower[chain_e]);
                ei.upper = std::min(ei.upper, _orig.upper[chain_e]);
                chain_costs.push_back(_orig.cost_function[chain_e]);
            }
            ei.cost_function = CostFunction::Sum(chain_costs.begin(), chain_costs.end());
            auto simp_e  = simp_.add_edge(ei);
            for (const auto chain_e: chain) {
                simp_edge_[chain_e] = simp_e;
            }
        }
    }
}


BiMDFResult BiMDF_Simplification::translate_solution(const BiMDFResult &_simp_result) const
{
    const auto &_simp_sol = *_simp_result.solution;
    auto orig_sol = std::make_unique<BiMDF::Solution>(orig_.g, 0);
    for (const auto e: orig_.g.edges()) {
        (*orig_sol)[e] = _simp_sol[simp_edge_[e]];
    }
    return {.solution = std::move(orig_sol),
            .cost = _simp_result.cost};
}

BiMDF_ConnectedComponents::BiMDF_ConnectedComponents(const BiMDF &_orig)
    : orig_(_orig)
    , node_cc_(_orig.g)
    , sub_node_(_orig.g)
    , sub_edge_(_orig.g)
{
    n_cc_ = lemon::connectedComponents(_orig.g, node_cc_);
    std::cout << "#CC: " << n_cc_ << std::endl;
    bimdfs_ = std::make_unique<BiMDF[]>(n_cc_);
    for (const auto n: _orig.g.nodes()) {
        auto cc = node_cc_[n];
        sub_node_[n] = bimdfs_[cc].add_node(_orig.demand[n]);
    }
    for (const auto e: _orig.g.edges()) {
        auto ei = _orig.get_edge_info(e);
        auto u_cc = node_cc_[ei.u];
        assert(u_cc == node_cc_[ei.v]);
        ei.u = sub_node_[ei.u];
        ei.v = sub_node_[ei.v];
        sub_edge_[e] =  bimdfs_[u_cc].add_edge(ei);

    }
    return;
}

BiMDFResult BiMDF_ConnectedComponents::translate_solutions(
            std::vector<BiMDFResult> const&_sols) const
{
    if (_sols.size() != n_cc_) {
        throw std::out_of_range("solution vector must be the same size as number of bimdfs!");
    }
    double total_cost = 0.;
    for (const auto &res: _sols) {
        total_cost += res.cost;
    }
    auto orig_sol = std::make_unique<BiMDF::Solution>(orig_.g, 0);
    for (const auto e: orig_.g.edges()) {
        auto u_cc = node_cc_[orig_.g.u(e)];
        assert(u_cc == node_cc_[orig_.g.v(e)]);
        const auto &sol = *_sols[u_cc].solution;
        (*orig_sol)[e] = sol[sub_edge_[e]];
    }

    return {.solution = std::move(orig_sol),
        .cost = total_cost};
}


} // namespace Satsuma
