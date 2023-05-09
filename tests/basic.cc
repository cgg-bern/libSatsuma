//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <libsatsuma/Problems/MCF.hh>

class BasicTest : public ::testing::Test { };

TEST_F(BasicTest, mcf_test)
{
    Satsuma::MCF mcf;
}

