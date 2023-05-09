//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "graphviz.hh"
#include "format.hh"

void save_matching_as_graphviz(
        std::ostream &s,
        Matching const& m)
{
    s << "graph {\n"
      << "rank = same;\n";
    const auto &g = m.g;
    for (const auto e: g.edges()) {
        std::string label = fmt_dbl(m.weight[e]);
            //<< ", u = " << mcf.upper[a]
        s << g.id(g.u(e))
          << " -- " << g.id(g.v(e))
          << " [label=\"" << label << "\"]\n";
    }
    s << "}\n";
}

void save_bmatching_as_graphviz(
        std::ostream &s,
        BMatching const& bm)
{
    s << "graph {\n"
      << "rank = same;\n";
    const auto &g = bm.g;

    for (const auto n: g.nodes()) {
        s << g.id(n) << "[label=\"" << g.id(n) << ": " << bm.degree[n] << "\"]\n";
    }

    for (const auto e: g.edges()) {
        std::string label = "(c=" + fmt_cap(bm.capacity[e]) + ", ω=" + fmt_dbl(bm.weight[e]) + ")";
            //<< ", u = " << mcf.upper[a]
        s << g.id(g.u(e))
          << " -- " << g.id(g.v(e))
          << " [label=\"" << label << "\"]\n";
    }
    s << "}\n";
}

void save_bimdf_as_graphviz(
        std::string filename,
        BiMDF const& bimdf)
{
    auto f = std::ofstream(filename);
    save_bimdf_as_graphviz(f, bimdf);
}
void save_bimcf_as_graphviz(
        std::string filename,
        BiMCF const& bimcf)
{
    auto f = std::ofstream(filename);
    save_bimcf_as_graphviz(f, bimcf);
}

void save_mcf_as_graphviz(
        std::string filename,
        MCF const& mcf)
{
    auto f = std::ofstream(filename);
    save_mcf_as_graphviz(f, mcf);
}

void save_bmatching_as_graphviz(
        std::string filename,
        BMatching const& bm)
{
    auto f = std::ofstream(filename);
    save_bmatching_as_graphviz(f, bm);
}

void save_matching_as_graphviz(
        std::string filename,
        Matching const& m)
{
    auto f = std::ofstream(filename);
    save_matching_as_graphviz(f, m);
}



void save_bimdf_as_graphviz(
        std::ostream &s,
        BiMDF const& bimdf)
{
    s << "digraph {\n"
      << "dir = LR;\n"
      << "rank = same;\n";
    const auto &g = bimdf.g;
    for (const auto e: g.edges()) {
        auto target = Satsuma::CostFunction::get_guess(bimdf.cost_function[e]);
        std::string label = "ℓ=" + fmt_dbl(target);
            //<< ", u = " << mcf.upper[e]
        s << g.id(g.u(e))
          << " -> " << g.id(g.v(e))
          << " [dir=both"
          << ",arrowhead=" << (bimdf.v_head[e] ? "\"normal\"" : "inv")
          << ",arrowtail=" << (bimdf.u_head[e] ? "\"normal\"" : "inv")
          << ",label=\"" << label << "\"]\n";
    }
    s << "}\n";
}
void save_bimcf_as_graphviz(
        std::ostream &s,
        BiMCF const& bimcf)
{
    s << "digraph {\n"
      << "dir = LR;\n"
      << "rank = same;\n";
    const auto &g = bimcf.g;
    for (const auto e: g.edges()) {
        std::string label = "(ω=" + fmt_dbl(bimcf.cost[e])
                          + ", u=" + fmt_cap(bimcf.upper[e]) + ")";
            //<< ", u = " << mcf.upper[e]
        s << g.id(g.u(e))
          << " -> " << g.id(g.v(e))
          << " [dir=both"
          << ",arrowhead=" << (bimcf.v_head[e] ? "\"normal\"" : "inv")
          << ",arrowtail=" << (bimcf.u_head[e] ? "\"normal\"" : "inv")
          << ",label=\"" << label << "\"]\n";
    }
    s << "}\n";
}

void save_mcf_as_graphviz(
        std::ostream &s,
        MCF const& mcf)
{
    s << "digraph {\n"
      << "rank = same;\n";
    const auto &g = mcf.g;
    for (const auto a: g.arcs()) {
        std::string label = "(ω=" + fmt_dbl(mcf.cost[a]/100.)
                          + ", u=" + fmt_cap(mcf.upper[a]) + ")";
            //<< ", u = " << mcf.upper[a]
        s << g.id(g.source(a))
          << " -> " << g.id(g.target(a))
          << " [label=\"" << label << "\"]\n";
    }
    s << "}\n";
}

