
#ifndef UTILS_UTILS_LOGGER_HPP_
#define UTILS_UTILS_LOGGER_HPP_

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
#include "utils_logger_impl.hpp"
#endif

namespace pressio{

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
using logger = spdlog::logger;
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
  logger->log(spdlog::level::info, "Initializing pressio logger");
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
  logger->log(spdlog::level::info, "Initializing pressio logger");
#endif
}

#endif//PRESSIO_ENABLE_TPL_PYBIND11


// Return an existing logger or nullptr if a logger with such name doesn't exist.
// example: spdlog::get("my_logger")->info("hello {}", "world");

template <class T= void>
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
std::shared_ptr<::pressio::logger> get(const std::string &name)
#else
void get(const std::string &name)
#endif
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  return spdlog::details::registry::instance().get(name);
#endif
}

template <class T= void>
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
std::shared_ptr<::pressio::logger> getLogger()
#else
void getLogger()
#endif
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  return get("pressioLogger");
#else
#endif
}

template <class T= void>
void finalize()
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
  getLogger()->log(spdlog::level::info, "Finalizing pressio logger");
  spdlog::shutdown();
#endif
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
#define PRESSIOLOGGER_ERROR(logger, ...) PRESSIO_LOGGER_NOSRCLOC_CALL(logger, spdlog::level::err, __VA_ARGS__)
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

#endif  // UTILS_UTILS_LOGGER_HPP_
