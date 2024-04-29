//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {
void write_bimdf(BiMDF const&_bimdf, std::string const &filename);
} // namespace Satsuma
