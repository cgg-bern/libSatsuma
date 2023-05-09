//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once


template<typename GraphT>
struct FigureGraph {
    explicit FigureGraph(const GraphT &g)
        : node_name{g}
        , node_label{g}
        , node_props{g}
    {}
    template<typename T> using NodeMap = typename GraphT::template NodeMap<T>;
    NodeMap<std::string> node_name, node_label;
    NodeMap<std::vector<std::string>> node_props;
};

template<typename GraphT = lemon::ListGraph>
struct FigureUndirGraph : public FigureGraph<GraphT>
{
    explicit FigureUndirGraph(const GraphT &g)
        : FigureGraph<GraphT>{g}
        , edge_label{g}
        , edge_props{g}
    {}
    template<typename T> using EdgeMap = typename GraphT::template EdgeMap<T>;
    EdgeMap<std::string> edge_label;
    EdgeMap<std::vector<std::string>> edge_props;
};

template<typename GraphT>
struct FigureDirGraph : public FigureGraph<GraphT>
{
    explicit FigureDirGraph(const GraphT &g)
        : FigureGraph<GraphT>{g}
        , arc_label{g}
        , arc_props{g}
    {}
    template<typename T> using ArcMap = typename GraphT::template ArcMap<T>;
    ArcMap<std::string> arc_label;
    ArcMap<std::vector<std::string>> arc_props;
};


