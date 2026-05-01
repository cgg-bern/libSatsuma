//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <libsatsuma/Problems/CostFunction.hh>

using namespace Satsuma;


class BasicTest : public ::testing::Test { };


TEST_F(BasicTest, find_optimal_integer) {
    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::QuadDeviation{.target = 5.7, .weight = 1.0}, 0, 20), 6);

    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::QuadDeviation{.target = 15.0, .weight = 1.0}, 0, 3), 3);

    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::QuadDeviation{.target = -5.0, .weight = 1.0}, 0, 20), 0);

    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::AbsDeviation{.target = 3.2, .weight = 1.0}, 0, 2), 2);

    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::AbsDeviation{.target = 7.9, .weight = 1.0}, 0, 20), 8);

    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::ScaleFactor{.target = 2.5, .weight = 1.0}, 0, 10), 3);

    struct MyFunc : CostFunction::BaseObjective {
        double operator()(double l) const override { return (l - 4.7)*(l - 4.7); }
        double get_guess() const override { return 4.7; }
    };
    EXPECT_EQ(CostFunction::find_optimal_integer(
        CostFunction::VirtualObjective{std::make_shared<MyFunc>()}, 0, 20), 5);

    {
        std::vector<CostFunction::Function> comps{
            CostFunction::QuadDeviation{.target = 3.0, .weight = 1.0},
            CostFunction::AbsDeviation{.target = 4.0, .weight = 1.0}
        };
        EXPECT_EQ(CostFunction::find_optimal_integer(
            CostFunction::Sum{comps.begin(), comps.end()}, 0, 20), 3);
    }

    {
        std::vector<CostFunction::Function> comps{
            CostFunction::QuadDeviation{.target = 2.0, .weight = 1.0},
            CostFunction::AbsDeviation{.target = 5.0, .weight = 1.0}
        };
        int r9 = CostFunction::find_optimal_integer(
            CostFunction::Sum{comps.begin(), comps.end()}, 0, 20);
        EXPECT_GE(r9, 1);
        EXPECT_LE(r9, 6);
    }
}
