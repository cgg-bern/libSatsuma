//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once
#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {

std::unique_ptr<BiMDF::Guess> make_guess(const BiMDF &bimdf);

} // namespace Satsuma
