//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include "tikz.hh"
#include "format.hh"

#include <map>

using std::to_string;

using GraphT = lemon::ListGraph;

//template<typename GraphT>
void save_nodes_as_tikz(
        std::ostream &s,
        GraphT const &g,
        FigureGraph<GraphT> const&fig)
{
    for (const auto n: g.nodes()) {
        const auto &name = fig.node_name[n];
        const auto &label = fig.node_label[n];

        s << "\t\\node[draw";
        for (const auto &prop: fig.node_props[n]) {
            s << "," << prop;
        }
        s << "] (" << name
          << ") at (p" << fig.node_name[n]
          << ") {$" << label << "$};\n";
    }
}

//template<typename GraphT>
std::unique_ptr<typename GraphT::template EdgeMap<int>>
compute_undirected_edge_multis(GraphT const &g)
{
    using Node = typename GraphT::Node;
    std::map<std::pair<Node, Node>, int> multi_edge_count, multi_edge_cur;
    for (const auto e: g.edges()) {
        ++multi_edge_count[std::make_pair(g.u(e), g.v(e))];
    }
    auto result = std::make_unique<typename GraphT::template EdgeMap<int>>(g);
    for (const auto e: g.edges())
    {
        auto u = g.u(e);
        auto v = g.v(e);
        int multicnt = multi_edge_count[std::make_pair(u, v)];
        int multicur = multi_edge_cur[std::make_pair(u, v)]++;
        //float mul = multicur - float(multicnt-1)/2;
        (*result)[e] = 1 + 2*multicur - multicnt;
    }
    return result;
}


//template<typename GraphT>
void save_graph_as_tikz(
        std::string filename,
        GraphT const &g,
        FigureUndirGraph<GraphT> const&fig,
        std::function<std::string(typename GraphT::Edge)> const &get_edge_type)
{
    auto s = std::ofstream(filename);
    //s << "\\tikzset{\n";

    save_nodes_as_tikz(s, g, fig);

    auto edge_multi = compute_undirected_edge_multis(g);

    for (const auto e: g.edges())
    {
        auto u = g.u(e);
        auto v = g.v(e);
        bool loop = u == v;

        const auto &label = fig.edge_label[e];

        s << "\t\\draw[" << get_edge_type(e);
        for (const auto &prop: fig.edge_props[e]) {
            s << "," << prop;
        }
        auto multi = (*edge_multi)[e];
        if (!loop && multi != 0) {
            s << ",multi=" + to_string(multi);
        }
        if (loop) {
            // https://tex.stackexchange.com/questions/222403/how-loop-edge-set-the-arrow-type
            s << ",every loop/.style={" << get_edge_type(e) << "}";
            //s << ",loop"; // TODO: due to complications with bidirected self-loops, better to call a tex macro instead of adding an attribute?
        }
        s <<  "] (" << fig.node_name[g.u(e)] <<") to";
        int loop_in=150 + multi*90;
        int loop_out=210 + multi*90;
        if (loop) {
            s << " [in=" << loop_in << ",out=" << loop_out << ",loop]";
        }
        if (label.size()) {
            s << " node[pos=.5,below,flowlabel] {$" << label << "$}";
        }
        s << "  (" << fig.node_name[g.v(e)] << ");\n";
    }

    // << "}\n";
}




void save_mcf_as_tikz(
        std::ostream &s,
        MCF const& mcf)
{
    const auto &g = mcf.g;

    s << "\\tikzset{\n";

    using Node = MCF::Node;
    auto is_plus = [&](Node n) {
        return (g.id(n) & 1) == 0;
    };

    std::stringstream plus, minus;

    for (const auto n: g.nodes()) {
        std::stringstream &stream = is_plus(n) ? plus : minus;
        stream << "\t\\node[flownode] (n" << g.id(n)
               << ") at (p" << g.id(n)
               << ") {$" << -mcf.supply[n] << "$};\n";
    }
    s << "pics/plus-nodes/.code={\n" << plus.str() << "},\n";
    s << "pics/minus-nodes/.code={\n" << minus.str() << "},\n";
    plus.str("");
    minus.str("");

    auto multi_index = [&](auto arc) {
        auto src = g.source(arc);
        auto dst = g.target(arc);
        return std::make_pair(std::min(src,dst), std::max(src,dst));
    };
    // count multi-edges (duplicated code, FIXME)
    std::map<std::pair<Node, Node>, int> multi_edge_count, multi_edge_cur;
    for (const auto a: g.arcs()) {
        ++multi_edge_count[multi_index(a)];
    }

    for (const auto a: g.arcs())
    {
        auto u = g.source(a);
        auto v = g.target(a);
        int multicnt = multi_edge_count[multi_index(a)];
        int multicur = multi_edge_cur[multi_index(a)]++;
        int mul = 1 + 2*multicur - multicnt;
        if (u < v) { mul = -mul;}
        auto cost = mcf.cost[a];
        bool free = cost == 0;
        auto dbl_cost = cost / double(1LL<<20);
        //std::cout << mcf.cost[a] << " - " << free << std::endl;
        auto upper = fmt_cap(mcf.upper[a], true);
        std::stringstream &stream = (is_plus(u) && is_plus(v)) ? plus : minus; // minus is in the back
        //std::string label = "c=" + fmt_dbl(mcf.cost[a]/100.)
        //                  + ", u=" + fmt_cap(mcf.upper[a]);
        //std::string label = fmt_cap(mcf.upper[a]) + " / "
        //                  + fmt_dbl(dbl_cost);
        stream << "\t\\draw[->,flowedge"
               <<  (free ? ",freearc" : "")
               <<  (multicnt > 1 ? ",multi=" + to_string(mul) : "")
               <<  (mcf.upper[a] == BiMCF::inf() ? ",infcap" : (",cap=" + upper ))
               <<  "] (n" << g.id(g.source(a)) <<") to";
        if (!free) {
            stream << " node[pos=.5,below,flowlabel] {$"<< fmt_dbl(dbl_cost) << "$}";
        }
        stream << "  (n" << g.id(g.target(a)) << ");\n";
    }
    s << "pics/plus-edges/.code={\n" << plus.str() << "},\n";
    s << "pics/minus-edges/.code={\n" << minus.str() << "},\n";
    s << "}\n";
}


#if 0
(
        std::string filename,
        GraphT const &g,
        FigureUndirGraph<GraphT> const&fig,
        std::function<std::string(typename GraphT::Edge)> const &get_edge_type)
#endif

