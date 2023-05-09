//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Problems/MCF.hh>

namespace Satsuma {


/// Solve via asymmetric Double cover.
/// Node demands must be even numbers.
/// If additionally, lower and upper bounds are even (or $\infty$), the solution is exact.
class BiMCF_to_MCF
{
public:
    enum class Method {
        HalfSymmetric,
        HalfAsymmetric,
        FullSymmetric,
        NotEven, // returns *double* of the half-integral solution!
        Default = HalfAsymmetric
    };
    struct Config {
        Method method = Method::Default;
        std::unique_ptr<MCF::NodeMap<BiMCF::Node>> *out_orig_node = nullptr;
        std::unique_ptr<MCF::NodeMap<bool>> *out_node_is_plus = nullptr;
    };
    BiMCF_to_MCF(BiMCF const &_bimcf,
                 Config const &_config);
    MCF const& mcf() const {return mcf_;}
    BiMCFResult translate_solution(const MCFResult &mcf_result) const;
private:
    BiMCF const& bimcf_;
    Method method_;
    MCF mcf_;
    MCF::ArcMap<BiMCF::Edge> orig_bimcf_edge_{mcf_.g};
    double costmul_ = 100;
};

} // namespace Satsuma
