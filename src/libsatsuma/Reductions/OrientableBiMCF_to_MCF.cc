
#include <libsatsuma/Reductions/OrientableBiMCF_to_MCF.hh>
#include <lemon/adaptors.h>
#include <cmath>
#include <cassert>

namespace Satsuma {

OrientedBiMCF::OrientedBiMCF(BiMCF const &_bimcf, Orientation const &ori)
    : bimcf_(_bimcf)
{
    size_t n_nodes = bimcf_.g.maxNodeId() + 1;
    size_t n_arcs = bimcf_.g.maxEdgeId() + 1;
    mcf_.g.reserveNode(n_nodes);
    mcf_.g.reserveArc(n_arcs);
    for (size_t i = 0; i < n_nodes; ++i) {
        auto n = mcf_.g.addNode();
        assert(mcf_.g.id(n) == i);
        mcf_.supply[n] = - bimcf_.demand[bimcf_.g.nodeFromId(i)];
    }
    for (size_t i = 0; i < n_arcs; ++i) {
        auto e = bimcf_.g.edgeFromId(i);
        assert(bimcf_.g.valid(e));
        bool u_head = bimcf_.u_head[e] ^ ori[bimcf_.g.u(e)];
        bool v_head = bimcf_.v_head[e] ^ ori[bimcf_.g.v(e)];
        if (!(u_head ^ v_head)) {
            std::cout
                << "uh " << bimcf_.u_head[e]
                << ", vh " << bimcf_.v_head[e]
                << ", ori u " << ori[bimcf_.g.u(e)]
                << ", ori v " << ori[bimcf_.g.v(e)]
                << std::endl;
            throw std::runtime_error("Invalid Orientation for arc " + std::to_string(i));
        }
        auto src = v_head ? bimcf_.g.u(e) : bimcf_.g.v(e);
        auto dst = u_head ? bimcf_.g.u(e) : bimcf_.g.v(e);
        auto src_id = bimcf_.g.id(src);
        auto dst_id = bimcf_.g.id(dst);

        auto new_e = mcf_.g.addArc(mcf_.g.nodeFromId(src_id),
                                   mcf_.g.nodeFromId(dst_id));
        orig_bimcf_edge_[new_e] = e;
        assert(mcf_.g.id(new_e) == i);
        mcf_.upper[new_e] = bimcf_.upper[e];
        mcf_.lower[new_e] = bimcf_.lower[e];
        mcf_.cost[new_e] = std::llround(bimcf_.cost[e]*costmul_);
    }
}
BiMCFResult OrientedBiMCF::translate_solution(const MCFResult &mcf_result) const
{
    auto bimcf_solp = std::make_unique<BiMCF::Solution>(bimcf_.g, 0);
    auto &bimcf_sol = *bimcf_solp;
    const auto &mcf_sol = *mcf_result.solution;
    size_t n_arcs = bimcf_.g.maxEdgeId() + 1;
    BiMCF::FlowScalar max_flow = 0;
    for (const auto arc: mcf_.g.arcs()) {
        auto bimcf_e = orig_bimcf_edge_[arc];
        auto val = mcf_sol[arc];
        bimcf_sol[bimcf_e] = val;
        max_flow = std::max(max_flow, val);
    }
    return {.solution = std::move(bimcf_solp),
            .cost = mcf_result.cost / costmul_,
            .max_flow = max_flow};

}

} // namespace Satsuma
