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

#ifdef PRESSIO_ENABLE_TPL_MPI
#include <mpi.h>
#endif
#include "utils/src/logger/spdlog/common.hpp"
#include "utils/src/logger/spdlog/spdlog.hpp"
#include "utils/src/logger/spdlog/sinks/stdout_color_sinks.hpp"
#include "utils/src/logger/spdlog/sinks/basic_file_sink.hpp"

#define PRESSIO_LOG_LEVEL_TRACE		0
#define PRESSIO_LOG_LEVEL_DEBUG		1
#define PRESSIO_LOG_LEVEL_INFO		2
#define PRESSIO_LOG_LEVEL_WARN		3
#define PRESSIO_LOG_LEVEL_ERROR		4
#define PRESSIO_LOG_LEVEL_CRITICAL	5
#define PRESSIO_LOG_LEVEL_OFF		6

#if defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_TRACE
#elif !defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_INFO
#endif

namespace pressio{

enum class logto{terminal, file, fileAndTerminal, terminalAndFile};

namespace log{

enum class level
{
  trace		= PRESSIO_LOG_LEVEL_TRACE,
  debug		= PRESSIO_LOG_LEVEL_DEBUG,
  info		= PRESSIO_LOG_LEVEL_INFO,
  warn		= PRESSIO_LOG_LEVEL_WARN,
  err		= PRESSIO_LOG_LEVEL_ERROR,
  critical	= PRESSIO_LOG_LEVEL_CRITICAL,
  off		= PRESSIO_LOG_LEVEL_OFF,
  n_levels
};

/////////////////////////////////////
/////////// begin IMPL //////////////
/////////////////////////////////////
namespace impl{

template<class T = void>
spdlog::level::level_enum pressioLogLevelToSpdlogLevel(::pressio::log::level en)
{
  switch(en)
    {
    case ::pressio::log::level::trace:
      return spdlog::level::level_enum::trace;
    case ::pressio::log::level::debug:
      return spdlog::level::level_enum::debug;
    case ::pressio::log::level::info:
      return spdlog::level::level_enum::info;
    case ::pressio::log::level::warn:
      return spdlog::level::level_enum::warn;
    case ::pressio::log::level::err:
      return spdlog::level::level_enum::err;
    case ::pressio::log::level::critical:
      return spdlog::level::level_enum::critical;
    case ::pressio::log::level::off:
      return spdlog::level::level_enum::off;
    case ::pressio::log::level::n_levels:
      return spdlog::level::level_enum::n_levels;
    default:
      throw std::runtime_error("Invalid pressio::log::level enum");
    };
}

template<class T = bool>
inline T mpiIsInitialized()
{
  bool value = false;
#if defined PRESSIO_ENABLE_TPL_MPI
  int flag = 0; MPI_Initialized( &flag );
  if (flag==1) value = true;
#endif
  return value;
}

template <typename T = int>
T myRank()
{
  T rank = 0;
#if defined PRESSIO_ENABLE_TPL_MPI
  if (mpiIsInitialized()){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
#endif
  return static_cast<T>(rank);
}

template <typename T = void>
std::string appendRankIfNeeded(std::string fin)
{
  if (mpiIsInitialized()){
    const auto rank = myRank();
    return fin + "_" + std::to_string(rank);
  }
  else
    return fin;
}

template <typename T = void>
std::shared_ptr<spdlog::logger> create(::pressio::logto en,
				       std::string fileName = "log.txt")
{
  const auto loggerName = "pressioLogger";
  const auto mpiRank = myRank();

  using tsink = spdlog::sinks::stdout_color_sink_mt;
  using fsink = spdlog::sinks::basic_file_sink_mt;

  if (en == ::pressio::logto::terminal)
  {
    auto terminal_sink = std::make_shared<tsink>();
    terminal_sink->set_level(spdlog::level::info);
    terminal_sink->setMpiRank(mpiRank);
    return std::make_shared<spdlog::logger>(loggerName,
					    spdlog::sinks_init_list({terminal_sink}));
  }
  else if (en == ::pressio::logto::fileAndTerminal or
	   en == ::pressio::logto::terminalAndFile)
  {
    if (fileName == "void"){
      throw std::runtime_error("Invalid filename to initialize logger.");
    }

    auto terminal_sink = std::make_shared<tsink>();
    auto file_sink     = std::make_shared<fsink>(appendRankIfNeeded(fileName), true);

    terminal_sink->setMpiRank(mpiRank);
    // default to info level both
    terminal_sink->set_level(spdlog::level::info);
    file_sink->set_level(spdlog::level::info);

    return std::make_shared<spdlog::logger>
      (loggerName, spdlog::sinks_init_list({file_sink, terminal_sink}));
  }
  else if (en == ::pressio::logto::file)
  {
    if (fileName == "void"){
      throw std::runtime_error("Invalid filename to initialize logger.");
    }

    auto file_sink = std::make_shared<fsink>(appendRankIfNeeded(fileName), true);
    file_sink->set_level(spdlog::level::info);
    return std::make_shared<spdlog::logger>(loggerName,
					    spdlog::sinks_init_list({file_sink}));
  }
  else{
    throw std::runtime_error("Invalid logto::enum value for logger");
  }
}

}//end namespace impl
/////////////////////////////////////
/////////// end IMPL ////////////////
/////////////////////////////////////

template<typename ...Args>
void initialize(::pressio::logto en, Args && ... args)
{
  auto logger = ::pressio::log::impl::create(en, std::forward<Args>(args)...);

  // // our pattern
  // logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%#] [fnc=%!] : %v");

  // the logger is defauled to min of trace, but not that the sinks
  // can have different levels. By using trace for the main logger
  // we make sure it is up to the sinks or the global define to set the min
  // level of output.
  logger->set_level(spdlog::level::trace);

  // set the singleton
  spdlog::set_default_logger(logger);
}

// Return an existing logger or nullptr if a logger with such name doesn't exist.
// example: spdlog::get("my_logger")->info("hello {}", "world");
template <class T= void>
std::shared_ptr<spdlog::logger> get(const std::string &name)
{
  return spdlog::details::registry::instance().get(name);
}

template <typename T = void>
void setVerbosity(std::initializer_list<::pressio::log::level> levels)
{
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
}

// Set global format string.
// example: spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %l : %v");
template <typename ...Args>
void setPattern(Args && ... args){
  spdlog::set_pattern(std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
inline void trace(const FormatString &fmt, Args&&...args){
  spdlog::default_logger()->trace(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
inline void debug(const FormatString &fmt, Args&&...args){
  spdlog::default_logger()->debug(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
inline void info(const FormatString &fmt, Args&&...args){
  spdlog::default_logger()->info(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
inline void warn(const FormatString &fmt, Args&&...args){
  spdlog::default_logger()->warn(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
inline void critical(const FormatString &fmt, Args&&...args){
  spdlog::default_logger()->critical(fmt, std::forward<Args>(args)...);
}

}} // namespace pressio::log


#define PRESSIO_LOGGER_CALL(logger, level, ...) \
  if (logger) (logger)->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, level, __VA_ARGS__)

///// TRACE /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_TRACE
#define PRESSIOLOG_TRACE(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::trace, __VA_ARGS__)
#else
#define PRESSIOLOG_TRACE(...) (void)0
#endif

///// DEBUG /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_DEBUG
#define PRESSIOLOG_DEBUG(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::debug, __VA_ARGS__)
#else
#define PRESSIOLOG_DEBUG(...) (void)0
#endif

///// INFO /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_INFO
#define PRESSIOLOG_INFO(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::info, __VA_ARGS__)
#else
#define PRESSIOLOG_INFO(...) (void)0
#endif

///// WARN /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_WARN
#define PRESSIOLOG_WARN(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::warn, __VA_ARGS__)
#else
#define PRESSIOLOG_WARN(...) (void)0
#endif

///// ERROR /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_ERROR
#define PRESSIOLOG_ERROR(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::error, __VA_ARGS__)
#else
#define PRESSIOLOG_ERROR(...) (void)0
#endif

///// CRITICAL /////
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_CRITICAL
#define PRESSIOLOG_CRITICAL(...) PRESSIO_LOGGER_CALL(spdlog::default_logger(), spdlog::level::critical, __VA_ARGS__)
#else
#define PRESSIOLOG_CRITICAL(...) (void)0
#endif

#endif  // UTILS_LOGGER_UTILS_LOGGER_HPP_
