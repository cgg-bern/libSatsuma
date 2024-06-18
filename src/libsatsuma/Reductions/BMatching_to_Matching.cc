//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Reductions/BMatching_to_Matching.hh>
#include <libsatsuma/Solvers/Matching.hh>
#include <limits>
#include <cassert>
#include <cmath>


namespace Satsuma {

const BMatching_to_Matching::Config BMatching_to_Matching::default_config = {};

BMatching_to_Matching::BMatching_to_Matching(const BMatching &_bmatching,
                                             Config const& _config)
    : bmatching_(_bmatching)
{
    if (_config.out_orig_node) {
        *_config.out_orig_node =  std::make_unique<Matching::NodeMap<BMatching::Node>>(matching_.g);
    }
    if (_config.out_node_num) {
        *_config.out_node_num =  std::make_unique<Matching::NodeMap<int>>(matching_.g);
    }
    if (_config.out_internode_edge) {
        *_config.out_internode_edge =  std::make_unique<Matching::NodeMap<BMatching::Edge>>(matching_.g);
    }

    //BMatching::GraphT::EdgeMap<int> first_edge_id(bmatching_.g, -1);
    BMatching::GraphT::NodeMap<int> first_node_id(bmatching_.g, -1);

    // TODO: special case for degree===1 - maybe in adding inter-edges?
    size_t n_nodes = 0;
    size_t n_edges = 0;
    for (auto bmp_node: bmatching_.g.nodes()) {
        first_node_id[bmp_node] = n_nodes;
        n_nodes += bmatching_.degree[bmp_node];
    }
    for (auto bmp_edge: bmatching_.g.edges())
    {
        //first_edge_id[bmp_edge] = n_edges;
        auto demand_u = bmatching_.degree[bmatching_.g.u(bmp_edge)];
        auto demand_v = bmatching_.degree[bmatching_.g.v(bmp_edge)];

        auto cap = bmatching_.capacity[bmp_edge];
        if (cap < BMatching::inf()) {
            n_edges += cap * (demand_u + demand_v);
            n_nodes += 2;
        }
        else {
            n_edges += demand_u * demand_v;
        }
    }
#if 0
    if (verbosity > 1)
    {
        std::cerr << "b-matching |V| = " << bmatching_.g.maxNodeId()+1
                  << ", |E| = " << bmatching_.g.maxEdgeId()+1
                  << std::endl
                  << "->matching |V| = " << n_nodes
                  << ", |E| = " << n_edges
                  << std::endl;
    }
#endif
    if (n_nodes > std::numeric_limits<int>::max()
            || n_edges > std::numeric_limits<int>::max())
    {
        throw std::runtime_error("graph too big, switch first_*_id to size_t?");
    }
    matching_.g.reserveNode(n_nodes);
    matching_.g.reserveEdge(n_edges);

    for (auto bmp_node: bmatching_.g.nodes()) {
        auto degree = bmatching_.degree[bmp_node];
        //std::cout << "XXXX" << matching_.g.maxNodeId()+1 << " : " << first_node_id[bmp_node] << std::endl;
        assert(matching_.g.maxNodeId() +1 == first_node_id[bmp_node]);
        for (int i = 0; i < degree; ++i) {
            auto n = matching_.g.addNode();
            if (_config.out_orig_node) {
                (**_config.out_orig_node)[n] = bmp_node;
            }
            if (_config.out_node_num) {
                (**_config.out_node_num)[n] = i;
            }
        }
    }
    for (auto bmp_edge: bmatching_.g.edges()) {
        auto u = bmatching_.g.u(bmp_edge);
        auto v = bmatching_.g.v(bmp_edge);
        auto demand_u = bmatching_.degree[u];
        auto demand_v = bmatching_.degree[v];
        auto cap = bmatching_.capacity[bmp_edge];
        Matching::CostScalar weight = std::llround(costmul_ * bmatching_.weight[bmp_edge]);

#if 0
        std::cout << "\tconnecting "
            << bmatching_.g.id(u) << " (deg " << demand_u << ") - "
            << bmatching_.g.id(v) << " (deg " << demand_v << ")"
            << ", cap = " << cap
            << std::endl;
#endif
        if (cap == BMatching::inf()) {
            cap = std::max(demand_u, demand_v);
        }

        if (cap < std::min(demand_u, demand_v))
        {
            for (int i = 0; i < cap; ++i) {
                auto nleft = matching_.g.addNode();
                auto nright = matching_.g.addNode();
                if (_config.out_orig_node) {
                    (**_config.out_orig_node)[nleft] = lemon::INVALID;
                    (**_config.out_orig_node)[nright] = lemon::INVALID;
                }
                if (_config.out_node_num) {
                    (**_config.out_node_num)[nleft] = 2*i;
                    (**_config.out_node_num)[nright] = 2*i+1;
                }
                if (_config.out_internode_edge) {
                    (**_config.out_internode_edge)[nleft] = bmp_edge;
                    (**_config.out_internode_edge)[nright] = bmp_edge;
                }
#if 0
                std::cout << "\tLR " << i << ": " << matching_.g.id(nleft)
                          << ", "  << matching_.g.id(nright) << std::endl;
#endif
                auto inter = matching_.g.addEdge(nleft, nright);
                matching_.weight[inter] = 0;
                orig_edge[inter] = lemon::INVALID;

                for (int j = 0; j < demand_u; ++j) {
                    auto mu = matching_.g.nodeFromId(first_node_id[u] + j);
                    auto e = matching_.g.addEdge(mu, nleft);
                    matching_.weight[e] = weight;
                    orig_edge[e] = bmp_edge;

                }
                for (int j = 0; j < demand_v; ++j) {
                    auto mv = matching_.g.nodeFromId(first_node_id[v] + j);
                    auto e = matching_.g.addEdge(mv, nright);
                    matching_.weight[e] = 0;
                    orig_edge[e] = lemon::INVALID;
                }
            }
        } else {
            for (int i = 0; i < demand_u; ++i) {
                auto mu = matching_.g.nodeFromId(first_node_id[u] + i);
                assert(matching_.g.id(mu) <= matching_.g.maxNodeId());
                for (int j = 0; j < demand_v; ++j) {
                    auto mv = matching_.g.nodeFromId(first_node_id[v] + j);
                    if (mv == mu) // a self-loop can never be part of a matching.
                        continue;
                    assert(matching_.g.id(mv) <= matching_.g.maxNodeId());
                    auto e = matching_.g.addEdge(mu, mv);
                    matching_.weight[e] = weight;
                    orig_edge[e] = bmp_edge;
                }
            }
        }
    }


}

BMatchingResult BMatching_to_Matching::translate_solution(const MatchingResult &matching_result) const
{
    const auto &msol = *matching_result.solution;

    auto sol = std::make_unique<BMatching::Solution>(bmatching_.g, 0);
    for (auto e: matching_.g.edges()) {
        if (orig_edge[e] != lemon::INVALID) {
            (*sol)[orig_edge[e]] += msol[e];
        }
    }
    if (!bmatching_.is_valid(*sol)) {
        throw std::runtime_error("B-matching solution kaput");
    }
    return {.solution = std::move(sol), .weight = matching_result.weight / costmul_};
}

} // namespace Satsuma
