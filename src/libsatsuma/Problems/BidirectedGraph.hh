//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <lemon/list_graph.h>

namespace Satsuma {

struct BidirectedGraph
{
    using GraphT = lemon::ListGraph;
    using Node = typename GraphT::Node;
    using Edge = typename GraphT::Edge;
    using Arc = typename GraphT::Arc;

    template<typename T> using EdgeMap = typename GraphT::EdgeMap<T>;
    template<typename T> using ArcMap = typename GraphT::ArcMap<T>;
    template<typename T> using NodeMap = typename GraphT::NodeMap<T>;

    Node add_node() {
        ++n_nodes_;
        return g.addNode();
    }

    Edge add_edge(Node u, Node v) {
        ++n_edges_;
        return g.addEdge(u, v);
    }

    size_t n_nodes() const {return n_nodes_;}
    size_t n_edges() const {return n_edges_;}

public:
    /// Public member, but do not add/remove edges & vertices!
    /// TODO: is a const& getter sufficient?
    GraphT g;

    // for an bidirected edge (u,v) e: u_head[e]:
    //     true:  e has head at u (u <-? v)
    //     false: e has tail at u (u >-? v)
    EdgeMap<bool> u_head {g};
    EdgeMap<bool> v_head {g};

private:
    size_t n_nodes_ = 0;
    size_t n_edges_ = 0;
};

} // namespace Satsuma
