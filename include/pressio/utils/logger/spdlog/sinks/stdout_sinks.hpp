// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_SINKS_STDOUT_SINKS_HPP_
#define UTILS_LOGGER_SPDLOG_SINKS_STDOUT_SINKS_HPP_

#include "../details/console_globals.hpp"
#include "../details/synchronous_factory.hpp"
#include "./sink.hpp"
#include <cstdio>

#ifdef _WIN32
#include "../details/windows_include.hpp"
#endif

namespace spdlog { namespace sinks {

template<typename ConsoleMutex>
class stdout_sink_base : public sink
{
public:
  using mutex_t = typename ConsoleMutex::mutex_t;

  explicit stdout_sink_base(FILE *file);

  ~stdout_sink_base() override = default;

  stdout_sink_base(const stdout_sink_base &other) = delete;
  stdout_sink_base(stdout_sink_base &&other) = delete;

  stdout_sink_base &operator=(const stdout_sink_base &other) = delete;
  stdout_sink_base &operator=(stdout_sink_base &&other) = delete;

  void log(const details::log_msg &msg) override;
  void flush() override;
  void set_pattern(const std::string &pattern) override;
  void set_formatter(std::unique_ptr<spdlog::formatter> sink_formatter) override;

protected:
  mutex_t &mutex_;
  FILE *file_;
  std::unique_ptr<spdlog::formatter> formatter_;
#ifdef _WIN32
  HANDLE handle_;
#endif // WIN32
};

template<typename ConsoleMutex>
class stdout_sink : public stdout_sink_base<ConsoleMutex>
{
public:
  stdout_sink() : stdout_sink_base<ConsoleMutex>(stdout){}
};

template<typename ConsoleMutex>
class stderr_sink : public stdout_sink_base<ConsoleMutex>
{
public:
  stderr_sink() : : stdout_sink_base<ConsoleMutex>(stderr){}
};

using stdout_sink_mt = stdout_sink<details::console_mutex>;
using stdout_sink_st = stdout_sink<details::console_nullmutex>;

using stderr_sink_mt = stderr_sink<details::console_mutex>;
using stderr_sink_st = stderr_sink<details::console_nullmutex>;

} // namespace sinks

// factory methods
template<typename Factory = spdlog::synchronous_factory>
std::shared_ptr<logger> stdout_logger_mt(const std::string &logger_name);

template<typename Factory = spdlog::synchronous_factory>
std::shared_ptr<logger> stdout_logger_st(const std::string &logger_name);

template<typename Factory = spdlog::synchronous_factory>
std::shared_ptr<logger> stderr_logger_mt(const std::string &logger_name);

template<typename Factory = spdlog::synchronous_factory>
std::shared_ptr<logger> stderr_logger_st(const std::string &logger_name);

} // namespace spdlog

#ifdef SPDLOG_HEADER_ONLY
#include "stdout_sinks-inl.hpp"
#endif
#endif  // UTILS_LOGGER_SPDLOG_SINKS_STDOUT_SINKS_HPP_
