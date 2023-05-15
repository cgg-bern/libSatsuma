//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

namespace Satsuma {

enum class MatchingSolver {
    Lemon,
    BlossomV,
    Default = Lemon // always available
};
} // namespace Satsuma
