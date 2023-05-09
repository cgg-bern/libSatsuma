//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Problems/MCF.hh>
#include <libsatsuma/Problems/Matching.hh>
#include <libsatsuma/Problems/BMatching.hh>
#include <fstream>

using Satsuma::BiMDF;
using Satsuma::BiMCF;
using Satsuma::MCF;
using Satsuma::Matching;
using Satsuma::BMatching;

void save_bimdf_as_graphviz(std::ostream &s, BiMDF const& bimdf);
void save_bimcf_as_graphviz(std::ostream &s, BiMCF const& bimcf);
void save_mcf_as_graphviz(std::ostream &s, MCF const& mcf);
void save_bmatching_as_graphviz(std::ostream &s, BMatching const& bm);
void save_matching_as_graphviz(std::ostream &s, Matching const& m);

void save_bimdf_as_graphviz(std::string filename, BiMDF const& bimdf);
void save_bimcf_as_graphviz(std::string filename, BiMCF const& bimcf);
void save_mcf_as_graphviz(std::string filename, MCF const& mcf);
void save_bmatching_as_graphviz(std::string filename, BMatching const& bm);
void save_matching_as_graphviz(std::string filename, Matching const& m);



