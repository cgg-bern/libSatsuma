//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Problems/BMatching.hh>

namespace Satsuma {

enum class DeviationLimitKind {
    EdgeFlow,
    NodeThroughflow,
    Default = NodeThroughflow
};

class BiMCF_to_BMatching
{
public:
    struct Config {
        int max_deviation = 2;
        DeviationLimitKind deviation_limit = DeviationLimitKind::Default;
        std::unique_ptr<BMatching::NodeMap<BiMCF::Node>> *out_orig_node = nullptr;
        std::unique_ptr<BMatching::EdgeMap<BiMCF::Edge>> *out_orig_edge = nullptr; // value set to INVALID for inter-edges
        std::unique_ptr<BMatching::NodeMap<bool>> *out_is_in_node = nullptr;
    };
    BiMCF_to_BMatching(BiMCF const& _mcf, Config const& _config);
    BMatching const& bmatching() const {return bmatching_;}
    BiMCFResult translate_solution(const BMatchingResult &bmatching_result) const;

private:
    BiMCF const& bimcf_;
    BMatching bmatching_;
    BMatching::GraphT::EdgeMap<int> bm_edge_id = {bimcf_.g, -1}; // TODO use orig_edge instead and expose it for users
};

} // namespace Satsuma
