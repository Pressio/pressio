
#ifndef UTILS_LOGGER_UTILS_LOGGER_IMPL_HPP_
#define UTILS_LOGGER_UTILS_LOGGER_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <initializer_list>
#include "./spdlog/tweakme.hpp"
#include "./spdlog/details/null_mutex.hpp"
#include "./spdlog/common.hpp"
#include "./spdlog/details/os.hpp"
#include "./spdlog/details/log_msg.hpp"
#include "./spdlog/details/log_msg_buffer.hpp"
#include "./spdlog/details/circular_q.hpp"
#include "./spdlog/formatter.hpp"
#include "./spdlog/details/fmt_helper.hpp"
#include "./spdlog/pattern_formatter.hpp"
#include "./spdlog/details/backtracer.hpp"
#include "./spdlog/sinks/sink.hpp"
#include "./spdlog/logger.hpp"
#include "./spdlog/details/registry.hpp"
#include "./spdlog/details/synchronous_factory.hpp"
#include "./spdlog/spdlog.hpp"
#include "./spdlog/sinks/stdout_color_sinks.hpp"
#include "./spdlog/sinks/basic_file_sink.hpp"

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
  else if (en == ::pressio::logto::fileAndTerminal)
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
