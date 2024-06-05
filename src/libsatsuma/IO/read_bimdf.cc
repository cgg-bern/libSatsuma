//  SPDX-FileCopyrightText: 2024 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/IO/write_bimdf.hh>
#include <fstream>
#include <cassert>

namespace Satsuma {


static CostFunction::Function read_cost(std::istream &s)
{
    char cost_type;
    s >> cost_type;
    if (cost_type == '0') {
        return  CostFunction::Zero();
    } else if (cost_type == 'A') {
        double target, weight;
        s >> target >> weight;
        return  CostFunction::AbsDeviation({.target = target, .weight = weight});
    } else if (cost_type == 'Q') {
        double target, weight;
        s >> target >> weight;
        return  CostFunction::QuadDeviation({.target = target, .weight = weight});
    } else if (cost_type == 'S') {
        double target, weight, eps;
        s >> target >> weight >> eps;
        return  CostFunction::ScaleFactor({
                .target = target,
                .weight = weight,
                .eps = eps});
    } else if (cost_type == '+') {
        size_t count;
        double target, weight, eps;
        s >> count;
        std::vector<CostFunction::Function> cfs;
        for (size_t i = 0; i < count; ++i) {
            cfs.push_back(read_cost(s));
        }
        return  CostFunction::Sum(cfs.begin(), cfs.end());
    } else {
        throw std::runtime_error(std::string("Satsuma::read_bimdf: unknown cost type '") + cost_type + "'");
    }
}

std::unique_ptr<BiMDF> read_bimdf(std::string const &filename, bool verbose)
{
    std::fstream s(filename);
    if (!s.good()) {
        throw std::runtime_error("Satsuma::read_bimdf: cannot open file for reading.");
    }
    std::string comment;
    std::getline(s, comment);
    if (verbose) {
        std::cout << "Reading BiMDF, file comment " << comment << std::endl;
    }
    auto bimdf_p = std::make_unique<BiMDF>();
    auto &bimdf = *bimdf_p;
    auto &g = bimdf.g;
    s << std::hexfloat;
    size_t n_nodes, n_edges;
    s >> n_nodes;
    if (verbose) {
        std::cout << " n_nodes: " << n_nodes << std::endl;
    }
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
        int demand;
        s >> demand;
        auto n = bimdf.add_node(demand);
        assert(g.id(n) == node_id);
    }

    s >> n_edges;
    if (verbose) {
        std::cout << " n_edges: " << n_edges << std::endl;
    }
    for (size_t edge_id = 0; edge_id < n_edges; ++edge_id)
    {
        BiMDF::EdgeInfo ei;
        size_t node_u, node_v;
        s >> node_u >> node_v;
        if (node_u >= n_nodes || node_v >= n_nodes) {
            throw std::runtime_error("Satsuma::read_bimdf: invalid node id '"
                    + std::to_string(std::max(node_u, node_v)) + "'");
        }
        ei.u = g.nodeFromId(node_u);
        ei.v = g.nodeFromId(node_v);
        s >> ei.u_head >> ei.v_head;
        s >> ei.lower >> ei.upper;

        ei.cost_function = read_cost(s);
        auto e = bimdf.add_edge(ei);
        assert(g.id(e) == edge_id);
    }
    if (!s.good()) {
        throw std::runtime_error("Satsuma::read_bimdf: file went bad, maybe truncated?");
    }
    return bimdf_p;
}

} // namespace Satsuma
