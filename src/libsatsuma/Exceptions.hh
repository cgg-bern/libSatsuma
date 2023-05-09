//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#pragma once

#include <stdexcept>

namespace Satsuma {

/// Exception base class, should not be used directly by Satsuma code.
class Exception : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// Problem is infeasible
class InfeasibleError : public Exception {
    using Exception::Exception;
};

/// Problem is not bounded below
class UnboundedError : public Exception {
    using Exception::Exception;
};

/// Satsuma bug
class InternalError : public Exception {
    using Exception::Exception;
};

}
