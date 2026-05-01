//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Config/Export.hh>
#pragma once
#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {
struct SATSUMA_EXPORT BiMDFLowerBoundResult {
    BiMDF::CostScalar cost;
    BiMDF::FlowScalar max_deviation;
};

/// Compute a lower bound for the cost using double cover
SATSUMA_EXPORT BiMDFLowerBoundResult bimdf_lower_bound(const BiMDF &bimdf, int initial_maxdev = 5);

} // namespace Satsuma
