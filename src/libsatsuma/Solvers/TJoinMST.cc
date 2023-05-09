//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/TJoinMST.hh>
#include <lemon/kruskal.h>
#include <lemon/dfs.h>
#include <lemon/adaptors.h>
#include <cassert>


namespace Satsuma {


namespace detail {
template<typename GraphT, typename Callback>
class BacktrackVisitor : public lemon::DfsVisitor<GraphT>
{
public:
    BacktrackVisitor(Callback _cb)
        : cb_(std::forward<Callback>(_cb))
    {}
    void leave(typename GraphT::Node const &n) {
        last_left_ = n;
    }
    void backtrack(typename GraphT::Arc const& arc) const {
        cb_(arc, last_left_);
    }
private:
    typename GraphT::Node last_left_;
    Callback cb_;
};
template<typename GraphT, typename Callback>
BacktrackVisitor<GraphT, Callback> make_backtrack_visitor(GraphT const&, Callback _cb) {
    return BacktrackVisitor<GraphT, Callback>(std::forward<Callback>(_cb));
}
} // namespace detail




TJoinResult solve_tjoin_mst(TJoin const& tjoin)
{
    using GraphT = TJoin::GraphT;
    using Arc = GraphT::Arc;
    using Edge = GraphT::Edge;
    using Node = GraphT::Node;

    const auto& g = tjoin.g;

    TJoin::EdgeMap<bool> in_msf{g, false};
    lemon::kruskal(g, tjoin.cost, in_msf);

    auto subgraph = lemon::FilterEdges(g, in_msf);

    TJoin::NodeMap<bool> parity{g};
    for (const auto n: g.nodes()) {
        parity[n] = tjoin.t[n];
    }

    TJoinResult res;
    res.solution = std::make_unique<TJoin::Solution>(g, false);
    auto &sol = *res.solution;

    auto vis = detail::make_backtrack_visitor(g, [&](
                                              Arc const &g_arc,
                                              Node last_left)
    {
            auto src = g.source(g_arc);
            auto tgt = g.target(g_arc);

#if 0
            std::cerr << "backtrack " << g.id(g_arc)
                      << ", src = " << g.id(src)
                      << ", tgt = " << g.id(tgt)
                      << std::endl;
#endif

            if (src != last_left) {
                assert(tgt == last_left);
                std::swap(src, tgt);
            }
            if (parity[src]) {
                Edge e {g_arc};
                sol[e] = true;
                res.cost += tjoin.cost[e];
                parity[src] = !parity[src];
                parity[tgt]= !parity[tgt];
            }
    });

    auto dfs = lemon::DfsVisit<decltype(subgraph),decltype(vis)>(subgraph, vis);
    dfs.run();

    return res;
}


} // namespace Satsuma
