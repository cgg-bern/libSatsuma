#pragma once
#include <libsatsuma/Problems/BiMDF.hh>

namespace Satsuma {

std::unique_ptr<BiMDF::Guess> make_guess(const BiMDF &bimdf);

} // namespace Satsuma
