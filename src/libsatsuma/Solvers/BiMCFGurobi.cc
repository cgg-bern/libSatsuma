//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMCFGurobi.hh>
#include <gurobi_c++.h>

#include <cmath>

namespace Satsuma {

BiMCFResult solve_bimcf_gurobi(BiMCF const &bimdf, BiMCF::Guess *x0)
{
    const auto &g = bimdf.g;
    auto solp = std::make_unique<BiMCF::Solution>(g);
    auto &sol = *solp;
    const auto n_edges = g.maxEdgeId() + 1;
    //const auto n_nodes = g.maxNodeId() + 1;
    try {
        auto env = GRBEnv();
        auto model = GRBModel(env);

        auto edge_vars = model.addVars(n_edges, GRB_INTEGER);

        auto edge_var = [&](BiMCF::Edge e) -> GRBVar& {
            return edge_vars[g.id(e)];
        };

        // edge variable bounds
        for (const auto e: g.edges()) {
            auto &ev = edge_var(e);
            ev.set(GRB_DoubleAttr_LB, bimdf.lower[e]);
            if (bimdf.upper[e] != BiMCF::inf()) {
                ev.set(GRB_DoubleAttr_UB, bimdf.upper[e]);
            }
            ev.set(GRB_DoubleAttr_Obj, bimdf.cost[e]);
            if (x0) {
                ev.set(GRB_DoubleAttr_Start, (*x0)[e]);
            }
        }

        // conservation constraints
        for (const auto n: g.nodes()) {
            GRBLinExpr expr;
            for (const auto e: g.incEdges(n)) {
                auto head = g.u(e) == n ? bimdf.u_head[e] : bimdf.v_head[e];
                if (head) {
                    expr += edge_var(e);
                } else {
                    expr -= edge_var(e);
                }
            }
            model.addConstr(expr, GRB_EQUAL, bimdf.demand[n]);
        }

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        model.optimize();
        const auto status = model.get(GRB_IntAttr_Status);
        if (status != GRB_OPTIMAL) {
            throw std::runtime_error("Gurobi solver failed");
        }

        const auto obj = model.get(GRB_DoubleAttr_ObjVal);
        //std::cout << "gurobi obj: " << obj << std::endl;

        // save results
        for (const auto e: g.edges()) {
            auto &ev = edge_var(e);
            sol[e] = std::llround(ev.get(GRB_DoubleAttr_X));
        }

        return {.solution = std::move(solp), .cost = obj};

    } catch(GRBException const &e) {
        std::cerr << "GRB Error code = " << e.getErrorCode()
                  << ": " << e.getMessage()
                  << std::endl;
        throw e;
    }
}

} // namespace Satsuma
