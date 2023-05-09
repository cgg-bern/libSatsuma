//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <span>

namespace Satsuma {

/// Split BiMDF problem into connected components
/// For an efficient implementation, we would use dynamic subgraph adapters,
/// but some of the current code relies on accessing the complete graph.
/// TODO PERF: do use dynamic subgraph adapters
class BiMDF_ConnectedComponents
{
public:
    using Node = BiMDF::Node;
    using Edge = BiMDF::Edge;
    template<typename T> using NodeMap = BiMDF::NodeMap<T>;
    template<typename T> using EdgeMap = BiMDF::EdgeMap<T>;

    BiMDF_ConnectedComponents(BiMDF const &_orig);
    std::span<BiMDF> bimdfs() const {return {bimdfs_.get(), bimdfs_.get()+n_cc_};}
    BiMDFResult translate_solutions(std::vector<BiMDFResult> const&_sols) const;
private:
    BiMDF const &orig_;
    NodeMap<size_t> node_cc_;
    NodeMap<Node> sub_node_; // node id in the corresponding cc subgraph
    EdgeMap<Edge> sub_edge_; // node id in the corresponding cc subgraph
    size_t n_cc_ = 0;
    std::unique_ptr<BiMDF[]> bimdfs_;
};

/// Create simplified BiMDF problem by collapsing demand-0 nodes
/// that have exactly one outgoing and one incoming arc.
class BiMDF_Simplification
{
public:
    using Node = BiMDF::Node;
    using Edge = BiMDF::Edge;
    template<typename T> using NodeMap = BiMDF::NodeMap<T>;
    template<typename T> using EdgeMap = BiMDF::EdgeMap<T>;

    BiMDF_Simplification(BiMDF const &_orig);
    BiMDF const& bimdf() const {return simp_;}
    BiMDFResult translate_solution(const BiMDFResult &_simp_result) const;
private:
    BiMDF const &orig_;
    BiMDF simp_;
    EdgeMap<Edge> simp_edge_; // map on orig BiMDF
};


} // namespace Satsuma
