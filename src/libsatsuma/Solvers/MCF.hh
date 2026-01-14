//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Config/Export.hh>

#include <libsatsuma/Problems/MCF.hh>

namespace Satsuma {

SATSUMA_EXPORT MCFResult solve_mcf_via_lemon_netsimp(const MCF &mcf);


} // namespace Satsuma
