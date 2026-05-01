//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Config/Export.hh>
#pragma once
#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {

SATSUMA_EXPORT std::unique_ptr<BiMDF::Guess> make_guess(const BiMDF &bimdf);

} // namespace Satsuma
