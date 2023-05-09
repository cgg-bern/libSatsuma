//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once
#include <libsatsuma/Problems/BidirectedGraph.hh>
#include <limits>

// TODO: move somewhere else, Structures/?


namespace Satsuma {
struct BiFlowGraph : public BidirectedGraph
{
    using GraphT = BidirectedGraph::GraphT;
    using Node = BidirectedGraph::Node;
    using Edge = BidirectedGraph::Edge;
    using Arc = BidirectedGraph::Arc;

    using FlowScalar = int;

    // TODO: rename Solution & Guess to Flow?
    using Solution = EdgeMap<FlowScalar>; /// A valid, conservative flow
    using Guess = EdgeMap<FlowScalar>;    /// A valid flow

    template<typename T> using EdgeMap = typename BidirectedGraph::EdgeMap<T>;
    template<typename T> using NodeMap = typename BidirectedGraph::NodeMap<T>;

    static FlowScalar inf() {return std::numeric_limits<FlowScalar>::max();}

    NodeMap<FlowScalar> demand {g};
    EdgeMap<FlowScalar> upper {g}; ///          f <= uppper
    EdgeMap<FlowScalar> lower {g}; /// lower <= f

    /// are conservation constraints met?
    bool is_valid(Solution const&sol) const;
    inline Node add_node(int _demand = 0) {
        Node n = BidirectedGraph::add_node();
        demand[n] = _demand;
        return n;
    }
    /// Apply the flow `f`, save resulting net node flow contents to `out`
    void apply_flow(Guess const& f, NodeMap<FlowScalar> &out) const;
    std::unique_ptr<NodeMap<FlowScalar>> apply_flow(Guess const& f) const;
};

} // namespace Satsuma
