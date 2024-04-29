//  SPDX-FileCopyrightText: 2024 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/IO/write_bimdf.hh>
#include <fstream>

namespace Satsuma {

void write_bimdf(BiMDF const&_bimdf, std::string const &filename)
{
    std::ofstream s(filename);
    if (!s.good()) {
        throw std::runtime_error("Satsuma::write_bimdf: cannot open file for writing.");
    }
    s << std::hexfloat;
    s << "# libSatsuma BiMDF file\n";
    auto n_nodes = _bimdf.n_nodes();
    auto n_edges = _bimdf.n_edges();
    const auto &g = _bimdf.g;
    s << _bimdf.n_nodes() << "\n";
    for (size_t node_id = 0; node_id < n_nodes; ++node_id)
    {
        auto n = g.nodeFromId(node_id);
        s << _bimdf.demand[n] << std::endl;
    }
    s << _bimdf.n_edges() << "\n";
    for (size_t edge_id = 0; edge_id < n_edges; ++edge_id)
    {
        auto e = g.edgeFromId(edge_id);
        s << g.id(g.u(e))
          << " " << g.id(g.v(e))
          << " " << (_bimdf.u_head[e] ? "1" : "0")
          << " " << (_bimdf.v_head[e] ? "1" : "0")
          << " " << _bimdf.lower[e]
          << " " << _bimdf.upper[e]
          //<< " " << (_bimdf.upper[e] == _bimdf.inf() ? "inf" : std::to_string(_bimdf.upper[e]))
          << " ";
        const auto &cost_func = _bimdf.cost_function[e];
        if (std::holds_alternative<CostFunction::Zero>(cost_func)) {
            s << "0";
        } else if (auto cost_absdev = std::get_if<CostFunction::AbsDeviation>(&cost_func)) {
            s << "A " << cost_absdev->target << " " << cost_absdev->weight;
        } else if (auto cost_quaddev = std::get_if<CostFunction::QuadDeviation>(&cost_func)) {
            s << "Q " << cost_quaddev->target << " " << cost_quaddev->weight;
        } else if (std::holds_alternative<CostFunction::VirtualObjective>(cost_func)) {
            // TODO: if we ever need this, have user supply custom serializer
            throw std::runtime_error("Satsuma::write_bimdf: Saving virtual cost functions is not supported.");
        } else if (std::holds_alternative<CostFunction::Sum>(cost_func)) {
            throw std::runtime_error("Satsuma::write_bimdf: Saving 'Sum' functions is not supported.'Sum' ");
        }
        s << "\n";
    }
}

} // namespace Satsuma
