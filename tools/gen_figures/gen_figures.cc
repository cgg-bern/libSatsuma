//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/MCF.hh>
#include <libsatsuma/Solvers/EvenBiMDF.hh>
#include <libsatsuma/Solvers/BiMDFDoubleCover.hh>
#include <libsatsuma/Solvers/BiMDFRefinement.hh>
#include <libsatsuma/Extra/Highlevel.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_MCF.hh>
#include <libsatsuma/Reductions/BiMCF_to_BMatching.hh>
#include <libsatsuma/Reductions/BMatching_to_Matching.hh>
#include <libsatsuma/Solvers/MCF.hh>
#include <libsatsuma/Solvers/TJoinMST.hh>
#include <map>
#include <fstream>
#include <iomanip>
#include <limits>
#include <lemon/kruskal.h>

#include "format.hh"
#include "graphviz.hh"
#include "libsatsuma/Problems/CostFunction.hh"
#include "tikz.hh"
#include "figure.hh"

using Satsuma::BiMDF;
using Satsuma::BiMCF;
using Satsuma::BiMCF_to_MCF;
using Satsuma::TJoin;
using Satsuma::MCF;
using Satsuma::Matching;
using Satsuma::BMatching;
using Satsuma::DeviationLimitKind;

namespace CostFunction = Satsuma::CostFunction;

using std::to_string;

void prepare_bimcf_figure(BiMCF const &bimcf,
                          FigureUndirGraph<BiMCF::GraphT> &fig)
{
    for (const auto n: bimcf.g.nodes()) {
        auto &props = fig.node_props[n];
        props = {"flownode"};
        fig.node_label[n] = fmt_cap(bimcf.demand[n]);

    }
    for (const auto e: bimcf.g.edges()) {
        auto &props = fig.edge_props[e];
        props = {"flowedge"};
        if (bimcf.cost[e] == 0) {
            props.push_back("freearc");
        } else {
            fig.edge_label[e] = fmt_dbl(bimcf.cost[e]);
        }
        auto upper = bimcf.upper[e];
        if (upper == bimcf.inf()) {
            props.push_back("infcap");
        } else {
            props.push_back("cap=" + to_string(upper));
        }
    }
}

int main(int argc, char *argv[])
{
    BiMDF bimdf;
    std::map<std::string, BiMDF::Node> nodes;
    std::map<std::string, BiMDF::Edge> edges;

    auto x = bimdf.add_node(); // boundary
    nodes["x"] = x;
    auto a = bimdf.add_node();
    nodes["a"] = a;
    auto b = bimdf.add_node();
    nodes["b"] = b;
    auto c = bimdf.add_node();
    nodes["c"] = c;

    using Abs = CostFunction::AbsDeviation;
    edges["x_x"] = bimdf.add_edge({
            .u=x, .v=x,
            .u_head=false, .v_head=false,
            .cost_function = Abs{.target=0, .weight = 0},
            });
    edges["x_a"] = bimdf.add_edge({
            .u=x, .v=a,
            .u_head=true, .v_head=true,
            .cost_function = Abs{.target=4, .weight = 1},
            .lower=1}); // TODO: display "lower" in bimdf figure
    edges["a_b"] = bimdf.add_edge({
            .u=a, .v=b,
            .u_head=false, .v_head=false,
            .cost_function = Abs{.target=.7, .weight = 1},
            .lower=0});
    edges["a_c"] = bimdf.add_edge({
            .u=a, .v=c,
            .u_head=false, .v_head=false,
            .cost_function = Abs{.target=.4, .weight = 1},
            .lower=0});
    edges["b_c"] = bimdf.add_edge({
            .u=b, .v=c,
            .u_head=true, .v_head=true,
            .cost_function = Abs{.target=.2, .weight = 1},
            .lower=0});

    auto figure_bimdf = FigureUndirGraph<BiMDF::GraphT>(bimdf.g);
    for (const auto &[name, e]: edges) {
        auto target = CostFunction::get_guess(bimdf.cost_function[e]);
        figure_bimdf.edge_label[e] = fmt_dbl(target);
        figure_bimdf.edge_props[e] = {"flowedge"};
    }
    for (const auto &[name, n]: nodes) {
        figure_bimdf.node_label[n] = name;
        figure_bimdf.node_name[n] = name;
        figure_bimdf.node_props[n] = {"flownode"};
    }
    save_as_tikz("bimdf.tex", bimdf, figure_bimdf);

    // direct rounding:
    BiMDF::Guess rounded_target{bimdf.g};
    for (const auto &e: bimdf.g.edges()) {
        auto target = CostFunction::get_guess(bimdf.cost_function[e]);
        auto rounded = static_cast<int>(std::lround(target));
        rounded_target[e] = std::max(bimdf.lower[e], rounded);
    }

    auto rounding = Satsuma::TJoinBasedRounding(bimdf, 3);
    auto rounding_result = rounding.solve();
    //auto guess = Satsuma::guess_for_even_rhs(bimdf, 3, Satsuma::EveningMode::MST);
    

    ////////// BEGIN T-Join/MST figure {
    const auto& tjoin = rounding.tjoin();

    TJoin::EdgeMap<bool> in_msf{tjoin.g, false};
    lemon::kruskal(tjoin.g, tjoin.cost, in_msf);

    auto tjoin_sol = solve_tjoin_mst(tjoin);

    auto figure_tjoin = FigureUndirGraph<TJoin::GraphT>(tjoin.g);
    for (const auto n: tjoin.g.nodes()) {
        figure_tjoin.node_name[n] = figure_bimdf.node_name[n];
        figure_tjoin.node_label[n] = figure_bimdf.node_label[n];
        figure_tjoin.node_props[n] = {"flownode",
            tjoin.t[n] ? "tjoin=1" : "tjoin=0" };
    }
    for (const auto e: tjoin.g.edges()) {
        figure_tjoin.edge_label[e] = fmt_dbl(tjoin.cost[e]);
        figure_tjoin.edge_props[e] = {
            "flowedge",
            in_msf[e] ? "msf=1" : "msf=0",
            (*tjoin_sol.solution)[e] ? "tjoin-sol=1" : "tjoin-sol=0"
        };
    }
    save_ugraph_as_tikz("tjoin.tex", tjoin.g, figure_tjoin);

    ////////// END T-Join/MST figure }

    //auto red_bimcf_even = Satsuma::BiMDF_to_BiMCF(bimdf, *guess, 2, true, true);
    std::unique_ptr<BiMCF::NodeMap<BiMDF::Node>> bimcf_even_orig_node;

    auto red_bimcf_even = Satsuma::BiMDF_to_BiMCF(bimdf, {
            .guess = *rounding_result.guess,
            .max_deviation = 2,
            .last_arc_uncapacitated = true,
            .even = true,
            .consolidate = true,
            .out_orig_node = &bimcf_even_orig_node});
    const auto &bimcf_even = red_bimcf_even.bimcf();

    auto figure_bimcf_even = FigureUndirGraph<BiMCF::GraphT>(bimcf_even.g);
    for (const auto n: bimcf_even.g.nodes()) {
        auto orig = (*bimcf_even_orig_node)[n];
        figure_bimcf_even.node_name[n] = figure_bimdf.node_name[orig];
        figure_bimcf_even.node_label[n] = figure_bimdf.node_label[orig];
    }

    prepare_bimcf_figure(bimcf_even, figure_bimcf_even);
    save_as_tikz("bimcf_dc.tex", bimcf_even, figure_bimcf_even);

    std::unique_ptr<MCF::NodeMap<BiMCF::Node>> dc_orig_node;
    std::unique_ptr<MCF::NodeMap<bool>> dc_node_is_plus;

    auto red_dc = Satsuma::BiMCF_to_MCF(bimcf_even, {
                                            .method = BiMCF_to_MCF::Method::HalfAsymmetric,
                                            //.method = BiMCF_to_MCF::Method::NotEven,
                                            .out_orig_node = &dc_orig_node,
                                            .out_node_is_plus = &dc_node_is_plus
                                        });
    const auto &mcf_dc = red_dc.mcf();

    save_mcf_as_graphviz("mcf_dc.dot", mcf_dc);
    save_mcf_as_tikz("mcf_dc.tex", mcf_dc);

    auto figure_dc = FigureDirGraph<MCF::GraphT>(mcf_dc.g);
    for (const auto n: mcf_dc.g.nodes()) {
        auto orig = (*dc_orig_node)[n];
        auto plus = (*dc_node_is_plus)[n];
        figure_dc.node_name[n] = figure_bimdf.node_name[orig] + (plus ? "p" : "m");
        figure_dc.node_label[n] = figure_bimdf.node_label[orig] + (plus ? "^+" : "^-");
    }

    auto dc_sol = solve_mcf_via_lemon_netsimp(mcf_dc);
#if 0
    for (const auto e: mcf_dc.g.arcs()) {
        std::cout << " sol"
                  << " e " <<mcf_dc.g.id(e)
                  << ": " <<mcf_dc.g.id(mcf_dc.g.source(e))
                  << " -> " <<mcf_dc.g.id(mcf_dc.g.target(e))
                  << ": " << (*dc_sol.solution)[e]
                 << std::endl;
    }
#endif

    auto bimcf_sol_dc = red_dc.translate_solution(dc_sol);
    auto bimdf_sol_dc = red_bimcf_even.translate_solution(bimcf_sol_dc);
    std::cout << "DC cost: " << bimdf.cost(*bimdf_sol_dc.solution) << std::endl;
    std::cout << "----- finished double cover stuff -----" << std::endl;

    std::unique_ptr<BiMCF::NodeMap<BiMDF::Node>> bimcf_bm_orig_node;
    auto red_bm_bimcf = Satsuma::BiMDF_to_BiMCF(bimdf, {
            .guess = *bimdf_sol_dc.solution,
            .max_deviation = 2,
            .last_arc_uncapacitated = false,
            .even = false,
            .consolidate = true,
            .out_orig_node = &bimcf_bm_orig_node});
    const auto &bimcf_bm = red_bm_bimcf.bimcf();

    auto figure_bimcf_bm = FigureUndirGraph<BiMCF::GraphT>(bimcf_bm.g);
    for (const auto n: bimcf_bm.g.nodes()) {
        auto orig = (*bimcf_bm_orig_node)[n];
        figure_bimcf_bm.node_name[n] = figure_bimdf.node_name[orig];
        //figure_bimcf_bm.node_label[n] = figure_bimdf.node_label[orig];
    }
    prepare_bimcf_figure(bimcf_bm, figure_bimcf_bm);
    save_as_tikz("bimcf_bm.tex", bimcf_bm, figure_bimcf_bm);
    // restore labels (a,b,c), we have overwritten them with demands
    for (const auto n: bimcf_bm.g.nodes()) {
        auto orig = (*bimcf_bm_orig_node)[n];
        figure_bimcf_bm.node_label[n] = figure_bimdf.node_label[orig];
    }

    std::unique_ptr<BMatching::NodeMap<BiMCF::Node>> bm_orig_node;
    std::unique_ptr<BMatching::NodeMap<bool>> bm_is_in_node;
    std::unique_ptr<BMatching::EdgeMap<BiMCF::Edge>> bm_orig_edge;

    auto red_bm = Satsuma::BiMCF_to_BMatching(bimcf_bm, {
                                                  .max_deviation = 2,
                                                  .deviation_limit = DeviationLimitKind::NodeThroughflow,
                                                  .out_orig_node = &bm_orig_node,
                                                  .out_orig_edge = &bm_orig_edge,
                                                  .out_is_in_node = &bm_is_in_node});
    const auto &bm = red_bm.bmatching();
    auto figure_bm = FigureUndirGraph<BMatching::GraphT>(bm.g);
    for (const auto n: bm.g.nodes()) {
        auto orig = (*bm_orig_node)[n];
        auto is_in = (*bm_is_in_node)[n];
        figure_bm.node_name[n] = figure_bimcf_bm.node_name[orig] + (is_in ? "-i" : "-o");
        figure_bm.node_label[n] = figure_bimcf_bm.node_label[orig] + (is_in ? "^i" : "^o");
        figure_bm.node_props[n].push_back("flownode");
        if (!is_in) {
            figure_bm.node_props[n].push_back("node on layer=top");
        }
    }
    for (const auto e: bm.g.edges()) {
        if (!(*bm_is_in_node)[bm.g.u(e)] && !(*bm_is_in_node)[bm.g.v(e)]) {
            figure_bm.edge_props[e].push_back("on layer=top");
        }

        figure_bm.edge_props[e].push_back("flowedge");
        if (bm.weight[e] == 0) {
            figure_bm.edge_props[e].push_back("freearc");
        } else {
            figure_bm.edge_label[e] = fmt_dbl(bm.weight[e]);
        }
        // order matters: add "cap" after "freearc" to get bold cap-2 free arcs!
        figure_bm.edge_props[e].push_back("cap=" + to_string(bm.capacity[e]));
#if 0
        if (bm.g.u(e) == bm.g.v(e)) {
            figure_bm.edge_props[e].push_back("loop");
        }
#endif
    }

    save_ugraph_as_tikz("bm.tex", bm.g, figure_bm);


    std::unique_ptr<Matching::NodeMap<BMatching::Node>> m_orig_node;
    std::unique_ptr<Matching::NodeMap<int>> m_node_num;
    std::unique_ptr<Matching::NodeMap<BMatching::Edge>> m_internode_edge;

    auto red_matching = Satsuma::BMatching_to_Matching(red_bm.bmatching(), {
                                                           .out_orig_node = &m_orig_node,
                                                           .out_node_num = &m_node_num,
                                                           .out_internode_edge = &m_internode_edge,
                                                           });
    const auto &matching = red_matching.matching();

    auto figure_m = FigureUndirGraph<Matching::GraphT>(matching.g);
    for (const auto n: matching.g.nodes()) {
        auto orig = (*m_orig_node)[n];
        auto node_num = (*m_node_num)[n];
        if (orig != lemon::INVALID) {
            figure_m.node_props[n].push_back("flownode");
            figure_m.node_name[n] = figure_bm.node_name[orig] 
                + "-" + to_string(node_num);
            figure_m.node_label[n] = figure_bm.node_label[orig] + "_{" + to_string(node_num) + "}";
        } else {
            figure_m.node_props[n].push_back("internode");
            auto e = (*m_internode_edge)[n];
            assert(e != lemon::INVALID);
            auto ename = figure_bm.node_name[bm.g.u(e)] + "_" + figure_bm.node_name[bm.g.v(e)];
            //auto ename = to_string(std::min(uid, vid)) + "-" + to_string(std::max(uid,vid));
            figure_m.node_name[n] = "inter-" + ename
                + "_" + to_string(matching.g.id(e)) // ename is not unique for our multi-graphs
                + "_" + to_string(node_num/2)
                + (node_num & 1 ? "_b" : "_a");

            figure_m.node_label[n] = "";
        }
    }

    save_ugraph_as_tikz("matching.tex", matching.g, figure_m);


    auto sol_matching = Satsuma::solve_matching(red_matching.matching());

    auto sol_bmatching = red_matching.translate_solution(sol_matching);
    auto sol_bimcf = red_bm.translate_solution(sol_bmatching);
    auto sol_bimdf = red_bm_bimcf.translate_solution(sol_bimcf);
    std::cout << "BiMDF cost via BM: " << bimdf.cost(*sol_bimdf.solution) << std::endl;

    //auto red_bm = Satsuma::BiMCF_to_BMatching(bimcf_bm, 2);
    auto bimdf_sol1 = solve_bimdf_matching(bimdf);
    auto bimdf_sol2 = solve_bimdf_matching(bimdf);
    //std::cout << "sol cost " << bimdf_sol.cost << std::endl;

    std::cout << "graph sizes:"
              << "\n\tbi-mcf for bm: " << format_size(bimcf_bm.g)
              << "\n\tb-matching:    " << format_size(red_bm.bmatching().g)
              << "\n\tmatching:      "   << format_size(red_matching.matching().g)
              << std::endl;
    bool different = false;
    std::cout << "bimdf edges:\n";
    for (const auto &e: edges) {
        auto x_dc = (*bimdf_sol_dc.solution)[e.second];
        auto x_bm1 = (*bimdf_sol1.result.solution)[e.second];
        auto x_bm2 = (*bimdf_sol2.result.solution)[e.second];
        if (x_dc != x_bm2) {
            different = true;
        }
        std::cout << "  " << e.first << ": "
            << "ℓ = " << CostFunction::get_guess(bimdf.cost_function[e.second])
            << ", g = " << (*rounding_result.guess)[e.second]
            << ", x_dc = " << x_dc
            << ", x_bm1 = " << x_bm1
            << ", x_bm2 = " << x_bm2
            << std::endl;
    }
    if (!different) {
        std::cout << "DC and BM solution are identical." << std::endl;
        //return 1;
    }
    std::cout << "bimdf nodes:\n";
    //BiMDF::EdgeMap<BiMDF::FlowScalar> net(bimdf.g);
    auto net_round = bimdf.apply_flow(rounded_target);
    auto net_mst = bimdf.apply_flow(*rounding_result.guess);
    for (const auto &n: nodes) {
        std::cout << "  " << n.first << ": "
            << "net (round) = " << (*net_round)[n.second]
            << ", net (mst)   = " << (*net_mst)[n.second]
            << std::endl;
    }

#if 0
    std::cout << "bimcf edges:\n";
    for (const auto &e: edges) {
        std::cout << "  " << e.first << ": "
            << "ℓ = " << bimdf.target[e.second]
            << ", g = " << (*guess)[e.second]
            << ", x_dc = " << (*bimcf_sol_dc.solution)[e.second]
            << std::endl;
    }

    std::cout << "bimcf-even nodes:\n";
    for (const auto n: bimcf_even.g.nodes()) {
        std::cout << "  " << bimcf_even.g.id(n) << ": "
            << "demand = " << bimcf_even.demand[n]
            << std::endl;
    }


    std::cout << "mcf-dc nodes:\n";
    for (const auto n: mcf_dc.g.nodes()) {
        std::cout << "  " << mcf_dc.g.id(n) << ": "
            << "demand = " << mcf_dc.supply[n]
            << std::endl;
    }
    std::cout << "DC cost: " << dc_sol.cost << std::endl;
    std::cout << "bimcf-even edges:\n";
    for (const auto a: mcf_dc.g.arcs()) {
        std::cout << "  " << mcf_dc.g.id(a) << ": "
            << mcf_dc.g.id(mcf_dc.g.source(a))
            << " -> " << mcf_dc.g.id(mcf_dc.g.target(a))
            << ", cost " << mcf_dc.cost[a]
            << ", u = " << mcf_dc.upper[a]
            << ", sol = " << (*dc_sol.solution)[a]
            << std::endl;
    }
#endif


    save_bimdf_as_graphviz("bimdf.dot", bimdf);
    save_bimcf_as_graphviz("bimcf_dc.dot", bimcf_even);
    save_bimcf_as_graphviz("bimcf_bm.dot", bimcf_bm);
    save_bmatching_as_graphviz("bm.dot", red_bm.bmatching());
    save_matching_as_graphviz("matching.dot", red_matching.matching());
    


    return 0;
}
