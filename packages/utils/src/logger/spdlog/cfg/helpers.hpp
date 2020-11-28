// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_CFG_HELPERS_HPP_
#define UTILS_LOGGER_SPDLOG_CFG_HELPERS_HPP_

#include "../common.hpp"
#include <unordered_map>

namespace spdlog {
namespace cfg {
namespace helpers {
//
// Init levels from given string
//
// Examples:
//
// set global level to debug: "debug"
// turn off all logging except for logger1: "off,logger1=debug"
// turn off all logging except for logger1 and logger2: "off,logger1=debug,logger2=info"
//
void load_levels(const std::string &txt);
} // namespace helpers

} // namespace cfg
} // namespace spdlog

#ifdef SPDLOG_HEADER_ONLY
#include "helpers-inl.hpp"
#endif // SPDLOG_HEADER_ONLY
#endif  // UTILS_LOGGER_SPDLOG_CFG_HELPERS_HPP_
