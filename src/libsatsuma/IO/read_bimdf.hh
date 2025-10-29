//  SPDX-FileCopyrightText: 2024 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Config/Export.hh>

namespace Satsuma {

SATSUMA_EXPORT
std::unique_ptr<BiMDF> read_bimdf(std::string const &filename, bool verbose=false);

} // namespace Satsuma
