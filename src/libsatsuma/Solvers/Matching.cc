//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/Matching.hh>
#include <libsatsuma/Config/Blossom5.hh>
#include <lemon/matching.h>

#if SATSUMA_HAVE_BLOSSOM5
#  include <blossom5/PerfectMatching.h>
#endif

namespace Satsuma {

MatchingResult solve_matching_via_lemon(Matching const &mp)
{
    //lemon::TimeReport tr("MWPM-lemon: ");
    auto mwpm = lemon::MaxWeightedPerfectMatching(mp.g, mp.weight);
    auto success = mwpm.run();
    if (!success) {
        throw std::runtime_error("Matching problem is infeasible");
    }
    //std::cout << "mpwm cost: " << -mwpm.matchingWeight() << std::endl;
    auto sol = std::make_unique<Matching::Solution>(mp.g);
    for (auto e: mp.g.edges()) {
        (*sol)[e] = mwpm.matching(e);
    }
    return {.solution = std::move(sol), .weight = mwpm.matchingWeight()};

}
#if !SATSUMA_HAVE_BLOSSOM5
MatchingResult solve_matching_via_blossomV(Matching const &mp)
{
    throw std::runtime_error("Satsuma was built without blossom-V support.");
}
#else
MatchingResult solve_matching_via_blossomV(Matching const &mp)
{
    const auto &g = mp.g;
    auto n_nodes =g.maxNodeId()+1;
    auto n_edges = g.maxEdgeId()+1;
    BlossomV::PerfectMatching pm{n_nodes, n_edges};

#define DEBUG_BLOSSOMV 0
#if DEBUG_BLOSSOMV
    std::vector<int> costs, edges;
    costs.reserve(n_nodes);
    edges.reserve(2*n_edges);
#endif

    // TODO: Blossom allows updates to the input for warmstart!
    for (auto e: g.edges()) {
        int cost = -mp.weight[e];
        auto u = g.id(g.u(e));
        auto v = g.id(g.v(e));
        assert(u != v);
        pm.AddEdge(u, v, cost);
#if DEBUG_BLOSSOMV
        costs.push_back(cost);
        edges.push_back(u);
        edges.push_back(v);
#endif
    }
#if DEBUG_BLOSSOMV
    pm.Save("/tmp/matching.dimacs");
#endif
    pm.options.verbose = DEBUG_BLOSSOMV;
    pm.Solve();

#if DEBUG_BLOSSOMV
    int res = BlossomV::CheckPerfectMatchingOptimality(n_nodes, n_edges, edges.data(), costs.data(), &pm);
    printf("check optimality: res=%d (%s)\n", res, (res==0) ? "ok" : ((res==1) ? "error" : "fatal error"));
    double cost = BlossomV::ComputePerfectMatchingCost(n_nodes, n_edges, edges.data(), costs.data(), &pm);
    printf("cost = %.1f\n", cost);
#endif

    auto sol = std::make_unique<Matching::Solution>(mp.g);
    int bv_edge_id = 0;
    for (auto e: g.edges()) {
        (*sol)[e] = pm.GetSolution(bv_edge_id++);
    }
    //std::cout << "n_nodes " << g.maxNodeId()+1 << ", result count " << count << std::endl;
    auto total_weight = mp.cost(*sol);
    //std::cout << "blossom-V cost " << total_weight << std::endl;
    assert(mp.is_perfect(*sol));
    return {.solution = std::move(sol), .weight = total_weight};
}
#endif

MatchingResult solve_matching(const Matching &mp, MatchingSolver solver)
{
   switch (solver)
   {
   case MatchingSolver::BlossomV:
       return solve_matching_via_blossomV(mp);
   case MatchingSolver::Lemon:
       return solve_matching_via_lemon(mp);
   default:
       throw std::runtime_error("Unknown MatchingSolver");
   }
}


} // namespace Satsuma
