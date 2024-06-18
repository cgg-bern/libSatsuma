//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Solvers/BiMDFGurobi.hh>
#include <libsatsuma/Solvers/BiMCFGurobi.hh>
#include <libsatsuma/Solvers/BiMDFGuess.hh>
#include <libsatsuma/Reductions/BiMDF_to_BiMCF.hh>
#include <gurobi_c++.h>
#include <cmath>
#include <cassert>

namespace Satsuma {

BiMDFResult solve_bimdf_gurobi(BiMDF const &bimdf, bool int_targets, BiMDF::Guess *x0)
{
    const auto &g = bimdf.g;
    auto solp = std::make_unique<BiMDF::Solution>(g);
    auto &sol = *solp;
    const auto n_edges = g.maxEdgeId() + 1;
    //const auto n_nodes = g.maxNodeId() + 1;
    try {
        auto env = GRBEnv();
        auto model = GRBModel(env);

        auto edge_vars = model.addVars(n_edges, GRB_INTEGER);

        auto edge_var = [&](BiMDF::Edge e) -> GRBVar& {
            assert(g.id(e) < n_edges);
            return edge_vars[g.id(e)];
        };

        GRBQuadExpr obj; // TODO PERF: for linear cost functions, does this give us the q. solver?

        // edge variable bounds
        for (const auto e: g.edges()) {
            auto &ev = edge_var(e);
            if (x0) {
                ev.set(GRB_DoubleAttr_Start, (*x0)[e]);
            }

            ev.set(GRB_DoubleAttr_LB, bimdf.lower[e]);
            //model.addConstr(ev, GRB_GREATER_EQUAL, bimdf.lower[e]);
            //model.addConstr(ev, GRB_GREATER_EQUAL, 1);

            if (bimdf.upper[e] != BiMDF::inf()) {
                ev.set(GRB_DoubleAttr_UB, bimdf.upper[e]);
                //model.addConstr(ev, GRB_LESS_EQUAL, bimdf.upper[e]);
            }
            const auto &cost_func = bimdf.cost_function[e];
            if (std::holds_alternative<CostFunction::Zero>(cost_func)) {
                continue;
            } else if (auto f = std::get_if<CostFunction::AbsDeviation>(&cost_func)) {
                if (f->weight == 0.) {
                    continue;
                }
                auto slack = model.addVar(0, GRB_INFINITY, f->weight,
                                          int_targets ? GRB_INTEGER : GRB_CONTINUOUS);
                const auto target = f->target;
                if (int_targets) {
                    if (x0) {
                        slack.set(GRB_DoubleAttr_Start, std::fabs((*x0)[e] - std::llround(target)));
                    }
                    model.addConstr(slack, GRB_GREATER_EQUAL, ev - std::llround(target));
                    model.addConstr(slack, GRB_GREATER_EQUAL, std::llround(target) - ev);
                } else {
                    if (x0) {
                        slack.set(GRB_DoubleAttr_Start, std::fabs((*x0)[e] - target));
                    }
                    model.addConstr(slack, GRB_GREATER_EQUAL, ev - target);
                    model.addConstr(slack, GRB_GREATER_EQUAL, target - ev);
                }
            } else if (auto f = std::get_if<CostFunction::QuadDeviation>(&cost_func)) {
                if (f->weight == 0) {
                    continue;
                }

                GRBLinExpr dev = ev - f->target;
                obj += f->weight * dev * dev;
            }
            else {
              assert(false);
            }
        }

        // conservation constraints
        for (const auto n: g.nodes()) {
            GRBLinExpr expr;
            for (const auto e: g.incEdges(n)) {
                auto head = (g.u(e) == n) ? bimdf.u_head[e] : bimdf.v_head[e];
                if (head) {
                    expr += edge_var(e);
                } else {
                    expr -= edge_var(e);
                }
            }
            model.addConstr(expr, GRB_EQUAL, bimdf.demand[n]);
        }

        model.setObjective(obj);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        //model.write("/tmp/flow.mps");
        model.optimize();
        const auto status = model.get(GRB_IntAttr_Status);
        if (status != GRB_OPTIMAL) {
            throw std::runtime_error("Gurobi solver failed");
        }

        const auto obj_val = model.get(GRB_DoubleAttr_ObjVal);
        //std::cout << "gurobi obj: " << obj_val << std::endl;

        // save results
        for (const auto e: g.edges()) {
            auto &ev = edge_var(e);
            sol[e] = std::llround(ev.get(GRB_DoubleAttr_X));
        }

        return {.solution = std::move(solp)};

    } catch(GRBException const &e) {
        std::cerr << "GRB Error code = " << e.getErrorCode()
                  << ": " << e.getMessage()
                  << std::endl;
        throw e;
    }
}

BiMDFResult solve_bimdf_gurobi_mcf(const BiMDF &bimdf, BiMDF::Guess *x0)
{
    // TODO: x0 support!
    auto guessp = make_guess(bimdf);
    const BiMDF::Guess &guess = *guessp;
    //auto dc_x0 = x0 ? nullptr : approximate_bimdf_doublecover(bimdf).solution;
    //auto red = BiMDF_to_BiMCF(bimdf, x0 ? *x0 : *dc_x0 , BiMDF_to_BiMCF::Objective::SimpleAbs);
    //auto red = BiMDF_to_BiMCF(bimdf, x0 ? *x0 : *dc_x0 , BiMDF_to_BiMCF::Objective::ExactAbsRefinement, 5);
    auto red = BiMDF_to_BiMCF(bimdf, {
            .guess = guess,
            .last_arc_uncapacitated = true,
            .even = false});

    auto sol_mcf = solve_bimcf_gurobi(red.bimcf());
    return red.translate_solution(sol_mcf);
}


} // namespace Satsuma
