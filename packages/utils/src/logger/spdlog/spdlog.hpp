// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

// spdlog main header file.
// see example.cpp for usage example

#ifndef UTILS_LOGGER_SPDLOG_SPDLOG_HPP_
#define UTILS_LOGGER_SPDLOG_SPDLOG_HPP_

//#include "common.hpp"
//#include "details/registry.hpp"
//#include "logger.hpp"
//#include "details/synchronous_factory.hpp"
//#include "pattern_formatter.hpp"
//#include <chrono>

namespace spdlog {

using default_factory = synchronous_factory;

// // Create and register a logger with a templated sink type
// // The logger's level, formatter and flush level will be set according the
// // global settings.
// //
// // Example:
// //   spdlog::create<daily_file_sink_st>("logger_name", "dailylog_filename", 11, 59);
// template<typename Sink, typename... SinkArgs>
// std::shared_ptr<spdlog::logger> create(std::string logger_name, SinkArgs &&...sink_args)
// {
//     return default_factory::create<Sink>(std::move(logger_name), std::forward<SinkArgs>(sink_args)...);
// }

// API for using default logger (stdout_color_mt),
// e.g: spdlog::info("Message {}", 1);
//
// The default logger object can be accessed using the spdlog::default_logger():
// For example, to add another sink to it:
// spdlog::default_logger()->sinks().push_back(some_sink);
//
// The default logger can replaced using spdlog::set_default_logger(new_logger).
// For example, to replace it with a file logger.
//
// IMPORTANT:
// The default API is thread safe (for _mt loggers), but:
// set_default_logger() *should not* be used concurrently with the default API.
// e.g do not call set_default_logger() from one thread while calling spdlog::info() from another.

template<class T = void>
std::shared_ptr<spdlog::logger> default_logger()
{
  return details::registry::instance().default_logger();
}

template<class T = void>
spdlog::logger *default_logger_raw()
{
  return details::registry::instance().get_default_raw();
}

template<class T>
void set_default_logger(std::shared_ptr<T> default_logger)
{
  details::registry::instance().set_default_logger(std::move(default_logger));
}

// // Initialize and register a logger,
// // formatter and flush level will be set according the global settings.
// //
// // Useful for initializing manually created loggers with the global settings.
// //
// // Example:
// //   auto mylogger = std::make_shared<spdlog::logger>("mylogger", ...);
// //   spdlog::initialize_logger(mylogger);

// template<class T>
// void initialize_logger(std::shared_ptr<T> logger)
// {
//   details::registry::instance().initialize_logger(std::move(logger));
// }

// // Return an existing logger or nullptr if a logger with such name doesn't
// // exist.
// // example: spdlog::get("my_logger")->info("hello {}", "world");
// template<class T = void>
// std::shared_ptr<logger> get(const std::string &name)
// {
//   return details::registry::instance().get(name);
// }

// Set global formatter. Each sink in each logger will get a clone of this object
template<class T>
void set_formatter(std::unique_ptr<T> formatter)
{
  details::registry::instance().set_formatter(std::move(formatter));
}

// Set global format string.
// example: spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %l : %v");
template<class T= void>
void set_pattern(std::string pattern, pattern_time_type time_type = pattern_time_type::local)
{
  set_formatter(std::unique_ptr<spdlog::formatter>(new pattern_formatter(std::move(pattern), time_type)));
}

// // enable global backtrace support
// template<class T = void>
// void enable_backtrace(size_t n_messages)
// {
//   details::registry::instance().enable_backtrace(n_messages);
// }

// // disable global backtrace support
// template<class T = void>
// void disable_backtrace()
// {
//   details::registry::instance().disable_backtrace();
// }

// // call dump backtrace on default logger
// template<class T = void>
// void dump_backtrace()
// {
//   default_logger_raw()->dump_backtrace();
// }

// // Get global logging level
// template<class T = void>
// level::level_enum get_level()
// {
//   return default_logger_raw()->level();
// }

// // Set global logging level
// template<class T = void>
// void set_level(level::level_enum log_level)
// {
//   details::registry::instance().set_level(log_level);
// }

// // Determine whether the default logger should log messages with a certain level
// template<class T = void>
// bool should_log(level::level_enum lvl)
// {
//   return default_logger_raw()->should_log(lvl);
// }

// Set global flush level
template<class T = void>
void flush_on(level::level_enum log_level)
{
  details::registry::instance().flush_on(log_level);
}

// Start/Restart a periodic flusher thread
// Warning: Use only if all your loggers are thread safe!
template<class T = void>
void flush_every(std::chrono::seconds interval)
{
  details::registry::instance().flush_every(interval);
}

// // Set global error handler
// template<class T = void>
// void set_error_handler(void (*handler)(const std::string &msg))
// {
//   details::registry::instance().set_error_handler(handler);
// }

// // Register the given logger with the given name
// template<class T = void>
// void register_logger(std::shared_ptr<logger> logger)
// {
//   details::registry::instance().register_logger(std::move(logger));
// }

// // Apply a user defined function on all registered loggers
// // Example:spdlog::apply_all([&](std::shared_ptr<spdlog::logger> l) {l->flush();});
// template<class T>
// void apply_all(const std::function<void(std::shared_ptr<T>)> &fun)
// {
//   details::registry::instance().apply_all(fun);
// }

// // Drop the reference to the given logger
// template<class T = void>
// void drop(const std::string &name)
// {
//   details::registry::instance().drop(name);
// }

// // Drop all references from the registry
// template<class T = void>
// void drop_all()
// {
//   details::registry::instance().drop_all();
// }

// // stop any running threads started by spdlog and clean registry loggers
// template<class T = void>
// void shutdown()
// {
//   details::registry::instance().shutdown();
// }

// // Automatic registration of loggers when using spdlog::create() or spdlog::create_async
// template<class T = void>
// void set_automatic_registration(bool automatic_registration)
// {
//   details::registry::instance().set_automatic_registration(automatic_registration);
// }

// template<typename FormatString, typename... Args>
// void log(source_loc source, level::level_enum lvl, const FormatString &fmt, Args&&...args)
// {
//   default_logger_raw()->log(source, lvl, fmt, std::forward<Args>(args)...);
// }

// template<typename FormatString, typename... Args>
// void log(level::level_enum lvl, const FormatString &fmt, Args&&...args)
// {
//     default_logger_raw()->log(source_loc{}, lvl, fmt, std::forward<Args>(args)...);
// }

template<typename FormatString, typename... Args>
void trace(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->trace(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
void debug(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->debug(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
void info(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->info(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
void warn(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->warn(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
void error(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->error(fmt, std::forward<Args>(args)...);
}

template<typename FormatString, typename... Args>
void critical(const FormatString &fmt, Args&&...args)
{
    default_logger_raw()->critical(fmt, std::forward<Args>(args)...);
}

template<typename T>
void log(source_loc source, level::level_enum lvl, const T &msg)
{
    default_logger_raw()->log(source, lvl, msg);
}

template<typename T>
void log(level::level_enum lvl, const T &msg)
{
    default_logger_raw()->log(lvl, msg);
}

template<typename T>
void trace(const T &msg)
{
    default_logger_raw()->trace(msg);
}

template<typename T>
void debug(const T &msg)
{
    default_logger_raw()->debug(msg);
}

template<typename T>
void info(const T &msg)
{
    default_logger_raw()->info(msg);
}

template<typename T>
void warn(const T &msg)
{
    default_logger_raw()->warn(msg);
}

template<typename T>
void error(const T &msg)
{
    default_logger_raw()->error(msg);
}

template<typename T>
void critical(const T &msg)
{
    default_logger_raw()->critical(msg);
}

} // namespace spdlog

// //
// // enable/disable log calls at compile time according to global level.
// //
// // define SPDLOG_ACTIVE_LEVEL to one of those (before including spdlog.h):
// // SPDLOG_LEVEL_TRACE,
// // SPDLOG_LEVEL_DEBUG,
// // SPDLOG_LEVEL_INFO,
// // SPDLOG_LEVEL_WARN,
// // SPDLOG_LEVEL_ERROR,
// // SPDLOG_LEVEL_CRITICAL,
// // SPDLOG_LEVEL_OFF
// //

// #define SPDLOG_LOGGER_CALL(logger, level, ...) (logger)->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, level, __VA_ARGS__)

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
// #define SPDLOG_LOGGER_TRACE(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::trace, __VA_ARGS__)
// #define SPDLOG_TRACE(...) SPDLOG_LOGGER_TRACE(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_TRACE(logger, ...) (void)0
// #define SPDLOG_TRACE(...) (void)0
// #endif

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
// #define SPDLOG_LOGGER_DEBUG(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::debug, __VA_ARGS__)
// #define SPDLOG_DEBUG(...) SPDLOG_LOGGER_DEBUG(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_DEBUG(logger, ...) (void)0
// #define SPDLOG_DEBUG(...) (void)0
// #endif

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_INFO
// #define SPDLOG_LOGGER_INFO(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::info, __VA_ARGS__)
// #define SPDLOG_INFO(...) SPDLOG_LOGGER_INFO(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_INFO(logger, ...) (void)0
// #define SPDLOG_INFO(...) (void)0
// #endif

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_WARN
// #define SPDLOG_LOGGER_WARN(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::warn, __VA_ARGS__)
// #define SPDLOG_WARN(...) SPDLOG_LOGGER_WARN(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_WARN(logger, ...) (void)0
// #define SPDLOG_WARN(...) (void)0
// #endif

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_ERROR
// #define SPDLOG_LOGGER_ERROR(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::err, __VA_ARGS__)
// #define SPDLOG_ERROR(...) SPDLOG_LOGGER_ERROR(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_ERROR(logger, ...) (void)0
// #define SPDLOG_ERROR(...) (void)0
// #endif

// #if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_CRITICAL
// #define SPDLOG_LOGGER_CRITICAL(logger, ...) SPDLOG_LOGGER_CALL(logger, spdlog::level::critical, __VA_ARGS__)
// #define SPDLOG_CRITICAL(...) SPDLOG_LOGGER_CRITICAL(spdlog::default_logger_raw(), __VA_ARGS__)
// #else
// #define SPDLOG_LOGGER_CRITICAL(logger, ...) (void)0
// #define SPDLOG_CRITICAL(...) (void)0
// #endif

// // #ifdef SPDLOG_HEADER_ONLY
// // #include "spdlog-inl.hpp"
// // #endif

#endif  // UTILS_LOGGER_SPDLOG_SPDLOG_HPP_
