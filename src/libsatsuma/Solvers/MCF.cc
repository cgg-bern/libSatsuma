//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/MCF.hh>
#include <libsatsuma/Exceptions.hh>
#include <lemon/network_simplex.h>

namespace Satsuma {

MCFResult solve_mcf_via_lemon_netsimp(const MCF &mcf)
{
    lemon::NetworkSimplex<MCF::GraphT, MCF::FlowScalar, MCF::CostScalar> solver{mcf.g};
    using LemonSolver = decltype(solver);
    solver.costMap(mcf.cost);
    solver.supplyMap(mcf.supply);
    solver.upperMap(mcf.upper);
    solver.lowerMap(mcf.lower);
    auto res = solver.run();
    if (res == LemonSolver::INFEASIBLE) {
        throw InfeasibleError("netsimp: flow problem infeasible");
    } else if (res == LemonSolver::UNBOUNDED) {
        throw UnboundedError("netsimp: flow problem unbounded");
    } else if (res == LemonSolver::OPTIMAL) {
    } else {
        throw InternalError("netsimp: unknown solver return value " + std::to_string(res));
    }
    auto sol = std::make_unique<MCF::Solution>(mcf.g);

    for (const auto a: mcf.g.arcs()) {
        (*sol)[a] = solver.flow(a);
    }
    return {.solution = std::move(sol), .cost = solver.totalCost()};
}



} // namespace Satsuma
