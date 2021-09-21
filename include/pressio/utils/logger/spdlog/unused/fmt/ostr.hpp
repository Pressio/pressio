//
// Copyright(c) 2016 Gabi Melman.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)
//

#ifndef UTILS_LOGGER_SPDLOG_UNUSED_FMT_OSTR_HPP_
#define UTILS_LOGGER_SPDLOG_UNUSED_FMT_OSTR_HPP_
//
// include bundled or external copy of fmtlib's ostream support
//

// #if !defined(SPDLOG_FMT_EXTERNAL)
// #ifdef SPDLOG_HEADER_ONLY
// #ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
// #endif
// #endif
#include "./bundled/ostream.hpp"
// #else
// #include <fmt/ostream.h>
// #endif
#endif  // UTILS_LOGGER_SPDLOG_UNUSED_FMT_OSTR_HPP_
