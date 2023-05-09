//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

inline std::string format_size(auto const&g) {
#if 0
  return std::format("|V| = {}, |E| = {}",
      g.maxNodeId()+1,
      g.maxEdgeId()+1);
#endif
  return "|V| = " + std::to_string(g.maxNodeId()+1)
      +", |E| = " + std::to_string(g.maxEdgeId()+1);
}

inline std::string fmt_cap(int val, bool latex=false) {
  if (val == std::numeric_limits<int>::max()) {
    return latex ? "\\infty" : "∞";
  }
    return std::to_string(val);
}

inline std::string fmt_dbl(double d, int precision=3) {
    if (std::fabs(d) <= 1e-30) d = 0;
  std::stringstream s;
  s << std::setprecision(precision) << d;
  return s.str();
}
