//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Problems/BiMCF.hh>

namespace Satsuma {

class BiMDF_to_BiMCF
{
public:
    struct Config {
        const BiMDF::Guess &guess;// result = guess + mcf flow
        BiFlowGraph::FlowScalar max_deviation = 2;// currently only honored for refinement mode
        bool last_arc_uncapacitated = true;
        bool even = false;
        bool consolidate = true;
        std::unique_ptr<BiMCF::NodeMap<BiMDF::Node>> *out_orig_node = nullptr;
    };
    BiMDF_to_BiMCF(const BiMDF &_mdf, Config const& config);

    BiMCF const& bimcf() const { return bimcf_;}

    /// XXX double_guess for half-integral optimum
    BiMDFResult translate_solution(const BiMCFResult &bimcf_res, bool double_guess=false) const;

private:

    BiMDF const& bimdf_;
    BiMDF::Guess guess_;

    BiMCF bimcf_;

    BiMCF::EdgeMap<bool> is_forward_;
    BiMCF::EdgeMap<int> mdf_edge_id_;

};
} // namespace Satsuma
