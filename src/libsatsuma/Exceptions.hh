//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <libsatsuma/Config/Export.hh>
#include <stdexcept>

namespace Satsuma {

/// Exception base class, should not be used directly by Satsuma code.
class SATSUMA_EXPORT Exception : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// Problem is infeasible
class SATSUMA_EXPORT InfeasibleError : public Exception {
    using Exception::Exception;
};

/// Problem is not bounded below
class SATSUMA_EXPORT UnboundedError : public Exception {
    using Exception::Exception;
};

/// Satsuma bug
class SATSUMA_EXPORT InternalError : public Exception {
    using Exception::Exception;
};

}
