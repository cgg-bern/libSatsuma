//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <variant>
#include <memory>
#include <vector>
#include <iterator>


namespace Satsuma::CostFunction {

struct Zero {
    double operator()(double) const {return 0.;}
    double get_guess() const {return guess;}
    double guess = 0.;
};

/// f(x) = weight * |x - target|
struct AbsDeviation {
    double target;
    double weight;
    double operator()(double l ) const {return weight * std::fabs((l-target));}
    double get_guess() const {return target;}
};

/// f(x) = weight * (x - target)^2
struct QuadDeviation {
    double target;
    double weight;
    double operator()(double l ) const {return weight * (l-target)*(l-target);}
    double get_guess() const {return target;}
};
/// f(x) = weight * max((x+eps)/(target+eps), (target+eps)/(x+eps))
struct ScaleFactor {
    double target;
    double weight;
    double eps = 0.1;
    double operator()(double l ) const {
        auto adj_tgt = target+eps;
        auto adj_l = l+eps;
        return weight * std::max(adj_l/adj_tgt, adj_tgt/adj_l);
    }
    double get_guess() const {return target;}
};

/// For VirtualObjective
struct BaseObjective {
    virtual ~BaseObjective();
    /// Evaluate f(x).
    virtual double operator()(double) const = 0;
    /// Get the minimum parameter value or another decent initial value.
    virtual double get_guess() const = 0;
};

/// User-defined convex(!) cost functions
struct VirtualObjective {
    double operator()(double l) const {return (*obj_)(l);}
    double get_guess() const {return obj_->get_guess();}
    /// shared_ptr so we can copy-assign this.
    std::shared_ptr<BaseObjective> obj_;
};

struct Sum;
//using BasicFunction = std::variant<Zero, AbsDeviation, QuadDeviation, VirtualObjective>;
using Function = std::variant<Zero, AbsDeviation, QuadDeviation, ScaleFactor, VirtualObjective, Sum>;
double cost(Function const&f, double _l);
double get_guess(Function const&f);

/// Sum of other types of cost functions
struct Sum {
    Sum(auto _begin, auto _end) {
        components_.reserve(std::distance(_begin, _end));
        auto min_guess = std::numeric_limits<double>::infinity();
        auto max_guess = -std::numeric_limits<double>::infinity();
        for (auto it = _begin; it != _end; ++it) {
            if (std::holds_alternative<Zero>(*it)) {
                continue;
            }
            auto g = Satsuma::CostFunction::get_guess(*it);
            min_guess = std::min(min_guess, g);
            max_guess = std::max(max_guess, g);
            // ensure a flat list of components:
            if (auto sump = std::get_if<Sum>(&*it)) {
                std::copy(sump->begin(), sump->end(), std::back_inserter(components_));
            } else {
                components_.push_back(*it);
            }
        }
        // Find integer guess in the range of component guesses.
        // PERF: Runtime could be improved by a variant of binary search, relying on convexity;
        //       For now, assume the range is small enough that it does not matter.
        auto best_cost = std::numeric_limits<double>::infinity();
        auto last_cost = std::numeric_limits<double>::infinity();
        for (auto g = std::llround(std::floor(min_guess));
             g <= std::llround(std::ceil(max_guess));
             ++g)
        {
            auto c = (*this)(g);
            if (c > last_cost) {
                // We are past the minimum of the convex function
                break;
            }
            last_cost = c;
            if (c < best_cost) {
                guess_ = g;
                best_cost = c;
            }
        }
    }
    double operator()(double l) const;
    double get_guess() const {return guess_;};
    auto begin() const {return components_.cbegin();}
    auto end() const {return components_.end();}
    size_t size() const {return components_.size();}
private:
    std::vector<Function> components_;
    double guess_ = 0.;
};

/// obj should be some variant of functions, e.g. Function
inline double cost(Function const&f, double _l) {
    return std::visit([_l](const auto &o) -> double {return o(_l);}, f);
}

/// obj should be some variant of functions, e.g. Function or BasicFunction
inline double get_guess(Function const&f) {
    return std::visit([](const auto &o) -> double {return o.get_guess();}, f);
}

inline double Sum::operator()(double l) const {
    double c = 0.;
    for (const auto &obj: components_) {
        c += cost(obj, l);
    }
    return c;
}
} // namespace Objective
