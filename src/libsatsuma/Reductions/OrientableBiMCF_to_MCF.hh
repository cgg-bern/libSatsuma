#pragma once

#include <libsatsuma/Problems/BiMCF.hh>
#include <libsatsuma/Problems/MCF.hh>
#include <libsatsuma/Solvers/OrientBinet.hh>

namespace Satsuma {

class OrientedBiMCF {
public:
    OrientedBiMCF(BiMCF const &_bimcf, Orientation const &_ori);
    MCF const& mcf() const {return mcf_;}
    BiMCFResult translate_solution(const MCFResult &mcf_result) const;
private:
    BiMCF const& bimcf_;
    MCF mcf_;
    MCF::ArcMap<BiMCF::Edge> orig_bimcf_edge_{mcf_.g};
    const double costmul_ = 1000;
};

} // namespace Satsuma
