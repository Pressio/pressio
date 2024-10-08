
#ifndef UTILS_LOGGER_UTILS_LOGGER_HPP_
#define UTILS_LOGGER_UTILS_LOGGER_HPP_

#include "pressio/macros.hpp"

#if defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
#include "./fmt/fmt.hpp"
#endif

namespace pressio{

enum class logto{
  terminal
#if defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
  ,file
  ,fileAndTerminal
#endif
};

namespace log{
enum class level
{
  trace		  = PRESSIO_LOG_LEVEL_TRACE,
  debug		  = PRESSIO_LOG_LEVEL_DEBUG,
  info		  = PRESSIO_LOG_LEVEL_INFO,
  warn		  = PRESSIO_LOG_LEVEL_WARN,
  err		    = PRESSIO_LOG_LEVEL_ERROR,
  critical	= PRESSIO_LOG_LEVEL_CRITICAL,
  off		    = PRESSIO_LOG_LEVEL_OFF,
  n_levels
};
}// end namespace pressio::log
}// end namespace pressio

// this impl has to be included here because it needs
// visibility of the enums above
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
#include "utils_logger_impl.hpp"
#endif

namespace pressio{

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
using logger_t = spdlog::logger;
#endif

namespace log{

template<typename ...Args>
void initialize(::pressio::logto en, Args && ... args)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
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

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
template <class T= void>
std::shared_ptr<::pressio::logger_t> get(const std::string &name)
{
  return spdlog::details::registry::instance().get(name);
}
#endif

template <class T= void>
void finalize()
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
  auto logger = ::pressio::log::get("pressioLogger");
  logger->log(spdlog::level::info, "Finalizing pressio logger");
  spdlog::shutdown();
#endif
}

template <typename T= void>
void setVerbosity(std::initializer_list<::pressio::log::level> levels)
{
#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF && defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)
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

}} // namespace pressio::log

#if !defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_TRACE
#define PRESSIOLOG_TRACE(...) (void)0
#endif

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_DEBUG
#define PRESSIOLOG_DEBUG(...) (void)0
#endif

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_INFO
#define PRESSIOLOG_INFO(...) (void)0
#endif

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_WARN
#define PRESSIOLOG_WARN(...) (void)0
#endif

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_ERROR
#define PRESSIOLOG_ERROR(...) (void)0
#endif

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL <= PRESSIO_LOG_LEVEL_CRITICAL
#define PRESSIOLOG_CRITICAL(...) (void)0
#endif

#else // defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)

#if PRESSIO_LOG_ACTIVE_MIN_LEVEL != PRESSIO_LOG_LEVEL_OFF
#define PRESSIO_LOGGER_CALL(logger, level, ...) \
  if (logger) (logger)->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, level, __VA_ARGS__);

#define PRESSIO_LOGGER_NOSRCLOC_CALL(logger, level, ...) \
  if (logger) (logger)->log(level, __VA_ARGS__);
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

#endif // defined(PRESSIO_ENABLE_INTERNAL_SPDLOG)



#endif  // UTILS_LOGGER_UTILS_LOGGER_HPP_
