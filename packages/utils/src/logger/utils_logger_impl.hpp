/*
//@HEADER
// ************************************************************************
//
// utils_logger_impl.hpp
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

#ifndef UTILS_LOGGER_UTILS_LOGGER_IMPL_HPP_
#define UTILS_LOGGER_UTILS_LOGGER_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <initializer_list>
#include "utils/src/logger/spdlog/tweakme.hpp"
#include "utils/src/logger/spdlog/details/null_mutex.hpp"
#include "utils/src/logger/spdlog/common.hpp"
#include "utils/src/logger/spdlog/details/os.hpp"
#include "utils/src/logger/spdlog/details/log_msg.hpp"
#include "utils/src/logger/spdlog/details/log_msg_buffer.hpp"
#include "utils/src/logger/spdlog/details/circular_q.hpp"
#include "utils/src/logger/spdlog/formatter.hpp"
#include "utils/src/logger/spdlog/details/fmt_helper.hpp"
#include "utils/src/logger/spdlog/pattern_formatter.hpp"
#include "utils/src/logger/spdlog/details/backtracer.hpp"
#include "utils/src/logger/spdlog/sinks/sink.hpp"
#include "utils/src/logger/spdlog/logger.hpp"
#include "utils/src/logger/spdlog/details/registry.hpp"
#include "utils/src/logger/spdlog/details/synchronous_factory.hpp"
#include "utils/src/logger/spdlog/spdlog.hpp"
#include "utils/src/logger/spdlog/sinks/stdout_color_sinks.hpp"
#include "utils/src/logger/spdlog/sinks/basic_file_sink.hpp"

namespace pressio{ namespace log{ namespace impl{

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

}}}//end namespace pressio::log::impl

#endif  // UTILS_LOGGER_UTILS_LOGGER_IMPL_HPP_
