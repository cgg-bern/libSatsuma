//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Config/Export.hh>
#pragma once

#include <libsatsuma/Problems/TJoin.hh>

namespace Satsuma {

SATSUMA_EXPORT TJoinResult solve_tjoin_mst(TJoin const& tjoin);

} // namespace Satsuma
