//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Problems/MCF.hh>
#include <libsatsuma/Problems/Matching.hh>
#include <libsatsuma/Problems/BMatching.hh>
#include <fstream>
#include <functional>
#include <lemon/list_graph.h>
#include "figure.hh"

using Satsuma::BiMDF;
using Satsuma::BiMCF;
using Satsuma::MCF;
using Satsuma::Matching;
using Satsuma::BMatching;

void save_graph_as_tikz(
        std::string filename,
        lemon::ListGraph const &g,
        FigureUndirGraph<lemon::ListGraph> const&fig,
        std::function<std::string(typename lemon::ListGraph::Edge)> const &get_edge_type);


inline void save_as_tikz(
        std::string filename,
        Satsuma::BidirectedGraph const &bimdf,
        FigureUndirGraph<lemon::ListGraph> const&fig)
{
    save_graph_as_tikz(filename, bimdf.g, fig, 
        [&](auto e){
        auto uh = bimdf.u_head[e];
        auto vh = bimdf.v_head[e];
        if (uh) {
            return vh? "<->" : "<-";
        } else {
            return vh? "->" : ">-<";
        }});
}

inline void save_ugraph_as_tikz(
        std::string filename,
        lemon::ListGraph const &g,
        FigureUndirGraph<lemon::ListGraph> const&fig)
{
    save_graph_as_tikz(filename, g, fig, 
        [](auto e){return "-";});
}


void save_mcf_as_tikz( std::ostream &s, MCF const& mcf);

inline void save_mcf_as_tikz(
        std::string filename,
        MCF const& mcf)
{
    auto f = std::ofstream(filename);
    save_mcf_as_tikz(f, mcf);
}

