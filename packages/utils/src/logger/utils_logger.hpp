/*
//@HEADER
// ************************************************************************
//
// utils_logger.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef UTILS_LOGGER_UTILS_LOGGER_HPP_
#define UTILS_LOGGER_UTILS_LOGGER_HPP_

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
#include "utils_logger_impl.hpp"
#endif

namespace pressio{

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
using logger = spdlog::logger;
#else
// need to make a trivial logger type so that instrumentation
// around pressio code can be a noop without needing to add preproc dirs
using logger = ::pressio::utils::empty;
#endif

namespace log{

#if !defined PRESSIO_ENABLE_TPL_PYBIND11
template<typename ...Args>
void initialize(::pressio::logto en, Args && ... args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  auto logger = ::pressio::log::impl::create(en, std::forward<Args>(args)...);
  // logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%#] [fnc=%!] : %v");

  // note that the sinks can have different levels. By using trace for the
  // main logger we make sure it is up to the sinks or the global
  // define to set the minlevel of output.
  logger->set_level(spdlog::level::trace);

  // set the singleton
  spdlog::set_default_logger(logger);
#endif
}

#else //PRESSIO_ENABLE_TPL_PYBIND11

// need this because for pybind11 cannot use variadic directly
void initialize(::pressio::logto en, std::string fileName = "log")
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  auto logger = ::pressio::log::impl::create(en, fileName);
  // logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%#] [fnc=%!] : %v");

  // note that the sinks can have different levels. By using trace for the
  // main logger we make sure it is up to the sinks or the global
  // define to set the minlevel of output.
  logger->set_level(spdlog::level::trace);

  // set the singleton
  spdlog::set_default_logger(logger);
#endif
}
#endif

// Return an existing logger or nullptr if a logger with such name doesn't exist.
// example: spdlog::get("my_logger")->info("hello {}", "world");

template <class T= void>
std::shared_ptr<::pressio::logger> get(const std::string &name)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  return spdlog::details::registry::instance().get(name);
#else
  return nullptr;
#endif
}

template <class T= void>
std::shared_ptr<::pressio::logger> getLogger()
{
  return get("pressioLogger");
}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename T> void setVerbosity(T levels)
#else
template <typename T= void> void setVerbosity(std::initializer_list<::pressio::log::level> levels)
#endif
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  auto logger = ::pressio::log::get("pressioLogger");
  std::vector<spdlog::sink_ptr> & sinks = logger->sinks();
  // make sure the number of levels is same as sinks
  if (levels.size() != sinks.size())
    throw std::runtime_error("log: number of levels to set != number of sinks");

  auto l = levels.begin();
  std::for_each(sinks.begin(), sinks.end(),
		[&l](spdlog::sink_ptr & s){
		  auto spdlog_l = impl::pressioLogLevelToSpdlogLevel(*l++);
		  s->set_level(spdlog_l);
		});
#endif
}

// Set global format string.
// example: spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %l : %v");
template <typename ...Args>
void setPattern(Args && ... args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::set_pattern(std::forward<Args>(args)...);
#endif
}

template<typename FormatString, typename... Args>
inline void trace(const FormatString &fmt, Args&&...args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::default_logger_raw()->trace(fmt, std::forward<Args>(args)...);
#endif
}

template<typename FormatString, typename... Args>
inline void debug(const FormatString &fmt, Args&&...args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::default_logger_raw()->debug(fmt, std::forward<Args>(args)...);
#endif
}

template<typename FormatString, typename... Args>
inline void info(const FormatString &fmt, Args&&...args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::default_logger_raw()->info(fmt, std::forward<Args>(args)...);
#endif
}

template<typename FormatString, typename... Args>
inline void warn(const FormatString &fmt, Args&&...args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::default_logger_raw()->warn(fmt, std::forward<Args>(args)...);
#endif
}

template<typename FormatString, typename... Args>
inline void critical(const FormatString &fmt, Args&&...args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  spdlog::default_logger_raw()->critical(fmt, std::forward<Args>(args)...);
#endif
}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename... Args>
// inline void print4py(spdlog::source_loc loc, spdlog::level::level_enum lvl, Args&&...args)
// {
//   {
//   pybind11::scoped_ostream_redirect stream
//     (
//      std::cout,                               // std::ostream&
//      pybind11::module_::import("sys").attr("stdout") // Python output
//      );

//   spdlog::memory_buf_t buf;
//   fmt::format_to(buf, std::forward<Args>(args)...);
//   pybind11::print(buf.data());
//   fmt::print(buf.data());
//   pybind11::print("\n");
//   //fmt::print(std::string(loc.filename) + " " + a + "\n");
//   }
// }
// #endif

}} // namespace pressio::log


#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
#define PRESSIO_LOGGER_CALL(logger, level, ...) \
  if (logger) (logger)->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, level, __VA_ARGS__);

#define PRESSIO_LOGGER_NOSRCLOC_CALL(logger, level, ...) if (logger) (logger)->log(level, __VA_ARGS__);
#endif

///// TRACE /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_TRACE
#define PRESSIOLOGGER_TRACE(logger, ...) PRESSIO_LOGGER_CALL(logger, spdlog::level::trace, __VA_ARGS__)
#define PRESSIOLOG_TRACE(...) PRESSIOLOGGER_TRACE(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_TRACE(...) (void)0
#endif

///// DEBUG /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_DEBUG
#define PRESSIOLOGGER_DEBUG(logger, ...) PRESSIO_LOGGER_CALL(logger, spdlog::level::debug, __VA_ARGS__)
#define PRESSIOLOG_DEBUG(...) PRESSIOLOGGER_DEBUG(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_DEBUG(...) (void)0
#endif

///// INFO /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_INFO
#define PRESSIOLOGGER_INFO(logger, ...) PRESSIO_LOGGER_NOSRCLOC_CALL(logger, spdlog::level::info, __VA_ARGS__)
#define PRESSIOLOG_INFO(...) PRESSIOLOGGER_INFO(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_INFO(...) (void)0
#endif

///// WARN /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_WARN
#define PRESSIOLOGGER_WARN(logger, ...) PRESSIO_LOGGER_NOSRCLOC_CALL(logger, spdlog::level::warn, __VA_ARGS__)
#define PRESSIOLOG_WARN(...) PRESSIOLOGGER_WARN(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_WARN(...) (void)0
#endif

///// ERROR /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_ERROR
#define PRESSIOLOGGER_ERROR(logger, ...) PRESSIO_LOGGER_NOSRCLOC_CALL(logger, spdlog::level::error, __VA_ARGS__)
#define PRESSIOLOG_ERROR(...) PRESSIOLOGGER_ERROR(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_ERROR(...) (void)0
#endif

///// CRITICAL /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_CRITICAL
#define PRESSIOLOGGER_CRITICAL(logger, ...) PRESSIO_LOGGER_NOSRCLOC_CALL(logger, spdlog::level::critical, __VA_ARGS__)
#define PRESSIOLOG_CRITICAL(...) PRESSIOLOGGER_CRITICAL(spdlog::default_logger(), __VA_ARGS__)
#else
#define PRESSIOLOG_CRITICAL(...) (void)0
#endif

#endif  // UTILS_LOGGER_UTILS_LOGGER_HPP_
