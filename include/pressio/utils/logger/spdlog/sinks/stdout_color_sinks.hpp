// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_SINKS_STDOUT_COLOR_SINKS_HPP_
#define UTILS_LOGGER_SPDLOG_SINKS_STDOUT_COLOR_SINKS_HPP_

// #ifdef _WIN32
// #include "../sinks/wincolor_sink.hpp"
// #else
#include "../sinks/ansicolor_sink.hpp"
//#endif

//#include "../details/synchronous_factory.hpp"

namespace spdlog { namespace sinks {

// #ifdef _WIN32
// using stdout_color_sink_mt = wincolor_stdout_sink_mt;
// using stdout_color_sink_st = wincolor_stdout_sink_st;
// using stderr_color_sink_mt = wincolor_stderr_sink_mt;
// using stderr_color_sink_st = wincolor_stderr_sink_st;
// #else
using stdout_color_sink_mt = ansicolor_stdout_sink_mt;
using stdout_color_sink_st = ansicolor_stdout_sink_st;
using stderr_color_sink_mt = ansicolor_stderr_sink_mt;
using stderr_color_sink_st = ansicolor_stderr_sink_st;
//#endif
} // namespace sinks

// template<typename Factory = spdlog::synchronous_factory>
// std::shared_ptr<logger> stdout_color_mt(const std::string &logger_name, color_mode mode = color_mode::automatic)
// {
//   return Factory::template create<sinks::stdout_color_sink_mt>(logger_name, mode);
// }

// template<typename Factory = spdlog::synchronous_factory>
// std::shared_ptr<logger> stdout_color_st(const std::string &logger_name, color_mode mode = color_mode::automatic)
// {
//   return Factory::template create<sinks::stdout_color_sink_st>(logger_name, mode);
// }

// template<typename Factory = spdlog::synchronous_factory>
// std::shared_ptr<logger> stderr_color_mt(const std::string &logger_name, color_mode mode = color_mode::automatic)
// {
//   return Factory::template create<sinks::stderr_color_sink_mt>(logger_name, mode);
// }

// template<typename Factory = spdlog::synchronous_factory>
// std::shared_ptr<logger> stderr_color_st(const std::string &logger_name, color_mode mode = color_mode::automatic)
// {
//   return Factory::template create<sinks::stderr_color_sink_st>(logger_name, mode);
// }
} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "stdout_color_sinks-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_SINKS_STDOUT_COLOR_SINKS_HPP_
