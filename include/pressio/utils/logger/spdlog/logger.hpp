// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_LOGGER_HPP_
#define UTILS_LOGGER_SPDLOG_LOGGER_HPP_

// Thread safe logger (except for set_error_handler())
// Has name, log level, vector of std::shared sink pointers and formatter
// Upon each log write the logger:
// 1. Checks if its log level is enough to log the message and if yes:
// 2. Call the underlying sinks to do the job.
// 3. Each sink use its own private copy of a formatter to format the message
// and send to its destination.
//
// The use of private formatter per sink provides the opportunity to cache some
// formatted data, and support for different format per sink.

//#include "./common.hpp"
//#include "./details/log_msg.hpp"
//#include "./details/backtracer.hpp"
//#include "./sinks/sink.hpp"
//#include "./pattern_formatter.hpp"
//#include <cstdio>
//#include <vector>

// #ifdef SPDLOG_WCHAR_TO_UTF8_SUPPORT
// #include "./details/os.hpp"
// #endif

#ifndef SPDLOG_NO_EXCEPTIONS
#define SPDLOG_LOGGER_CATCH()				\
  catch (const std::exception &ex)			\
    {							\
      err_handler_(ex.what());				\
    }							\
  catch (...)						\
    {							\
      err_handler_("Unknown exception in logger");	\
    }
#else
#define SPDLOG_LOGGER_CATCH()
#endif

namespace spdlog {

class logger
{
public:
  // Empty logger
  explicit logger(std::string name)
    : name_(std::move(name))
    , sinks_()
  {}

  // Logger with range on sinks
  template<typename It>
  logger(std::string name, It begin, It end)
    : name_(std::move(name))
    , sinks_(begin, end)
  {}

  // Logger with single sink
  logger(std::string name, sink_ptr single_sink)
    : logger(std::move(name), {std::move(single_sink)})
  {}

  // Logger with sinks init list
  logger(std::string name, sinks_init_list sinks)
    : logger(std::move(name), sinks.begin(), sinks.end())
  {}

  virtual ~logger() = default;

  logger(const logger &other)
    : name_(other.name_)
    , sinks_(other.sinks_)
    , level_(other.level_.load(std::memory_order_relaxed))
    , flush_level_(other.flush_level_.load(std::memory_order_relaxed))
    , custom_err_handler_(other.custom_err_handler_)
    , tracer_(other.tracer_)
  {}

  logger(logger &&other) SPDLOG_NOEXCEPT
    : name_(std::move(other.name_)),
      sinks_(std::move(other.sinks_)),
      level_(other.level_.load(std::memory_order_relaxed)),
      flush_level_(other.flush_level_.load(std::memory_order_relaxed)),
      custom_err_handler_(std::move(other.custom_err_handler_)),
      tracer_(std::move(other.tracer_))
  {}

  logger &operator=(logger other) SPDLOG_NOEXCEPT
  {
    this->swap(other);
    return *this;
  }

public:
  void swap(spdlog::logger &other) SPDLOG_NOEXCEPT
  {
    name_.swap(other.name_);
    sinks_.swap(other.sinks_);

    // swap level_
    auto other_level = other.level_.load();
    auto my_level = level_.exchange(other_level);
    other.level_.store(my_level);

    // swap flush level_
    other_level = other.flush_level_.load();
    my_level = flush_level_.exchange(other_level);
    other.flush_level_.store(my_level);

    custom_err_handler_.swap(other.custom_err_handler_);
    std::swap(tracer_, other.tracer_);
  }

  // // FormatString is a type derived from fmt::compile_string
  // template<
  //   typename FormatString,
  //   typename std::enable_if<fmt::is_compile_string<FormatString>::value, int>::type = 0,
  //   typename... Args
  //   >
  // void log(source_loc loc,
  // 	   level::level_enum lvl,
  // 	   const FormatString &fmt,
  // 	   Args&&...args)
  // {
  //   log_(loc, lvl, fmt, std::forward<Args>(args)...);
  // }

  // FormatString is NOT a type derived from fmt::compile_string but
  // is a string_view_t or can be implicitly converted to one
  template<typename... Args>
  void log(source_loc loc, level::level_enum lvl, string_view_t fmt, Args&&...args)
  {
    log_(loc, lvl, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void log(level::level_enum lvl, const FormatString &fmt, Args&&...args)
  {
    log(source_loc{}, lvl, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void trace(const FormatString &fmt, Args&&...args)
  {
    log(level::trace, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void debug(const FormatString &fmt, Args&&...args)
  {
    log(level::debug, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void info(const FormatString &fmt, Args&&...args)
  {
    log(level::info, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void warn(const FormatString &fmt, Args&&...args)
  {
    log(level::warn, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void error(const FormatString &fmt, Args&&...args)
  {
    log(level::err, fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void critical(const FormatString &fmt, Args&&...args)
  {
    log(level::critical, fmt, std::forward<Args>(args)...);
  }

  // template<typename T>
  // void log(level::level_enum lvl, const T &msg)
  // {
  //   log(source_loc{}, lvl, msg);
  // }

  // // T can be statically converted to string_view and isn't a fmt::compile_string
  // template<
  //   class T,
  //   typename std::enable_if<
  //     std::is_convertible<const T &, spdlog::string_view_t>::value && !fmt::is_compile_string<T>::value, int>::type = 0>
  // void log(source_loc loc, level::level_enum lvl, const T &msg)
  // {
  //   log(loc, lvl, string_view_t{msg});
  // }

  // void log(log_clock::time_point log_time, source_loc loc, level::level_enum lvl, string_view_t msg)
  // {
  //   bool log_enabled = should_log(lvl);
  //   bool traceback_enabled = tracer_.enabled();
  //   if (!log_enabled && !traceback_enabled)
  //     {
  // 	return;
  //     }

  //   details::log_msg log_msg(log_time, loc, name_, lvl, msg);
  //   log_it_(log_msg, log_enabled, traceback_enabled);
  // }

  // void log(source_loc loc, level::level_enum lvl, string_view_t msg)
  // {
  //   bool log_enabled = should_log(lvl);
  //   bool traceback_enabled = tracer_.enabled();
  //   if (!log_enabled && !traceback_enabled)
  //     {
  // 	return;
  //     }

  //   details::log_msg log_msg(loc, name_, lvl, msg);
  //   log_it_(log_msg, log_enabled, traceback_enabled);
  // }

  // void log(level::level_enum lvl, string_view_t msg)
  // {
  //   log(source_loc{}, lvl, msg);
  // }

  // // T cannot be statically converted to string_view or wstring_view
  // template<
  //   class T,
  //   typename std::enable_if<!std::is_convertible<const T &, spdlog::string_view_t>::value &&
  // 			    !is_convertible_to_wstring_view<const T &>::value,
  // 			    int>::type = 0>
  // void log(source_loc loc, level::level_enum lvl, const T &msg)
  // {
  //   log(loc, lvl, "{}", msg);
  // }

  template<typename T>
  void trace(const T &msg)
  {
    log(level::trace, msg);
  }

  template<typename T>
  void debug(const T &msg)
  {
    log(level::debug, msg);
  }

  template<typename T>
  void info(const T &msg)
  {
    log(level::info, msg);
  }

  template<typename T>
  void warn(const T &msg)
  {
    log(level::warn, msg);
  }

  template<typename T>
  void error(const T &msg)
  {
    log(level::err, msg);
  }

  template<typename T>
  void critical(const T &msg)
  {
    log(level::critical, msg);
  }

// #ifdef SPDLOG_WCHAR_TO_UTF8_SUPPORT
// #ifndef _WIN32
// #error SPDLOG_WCHAR_TO_UTF8_SUPPORT only supported on windows
// #else

//   template<typename... Args>
//   void log(source_loc loc, level::level_enum lvl, wstring_view_t fmt, Args&&...args)
//   {
//     bool log_enabled = should_log(lvl);
//     bool traceback_enabled = tracer_.enabled();
//     if (!log_enabled && !traceback_enabled)
//       {
// 	return;
//       }
//     SPDLOG_TRY
//       {
// 	// format to wmemory_buffer and convert to utf8
// 	fmt::wmemory_buffer wbuf;
// 	fmt::format_to(wbuf, fmt, std::forward<Args>(args)...);

// 	memory_buf_t buf;
// 	details::os::wstr_to_utf8buf(wstring_view_t(wbuf.data(), wbuf.size()), buf);
// 	details::log_msg log_msg(loc, name_, lvl, string_view_t(buf.data(), buf.size()));
// 	log_it_(log_msg, log_enabled, traceback_enabled, mpiRank_);
//       }
//     SPDLOG_LOGGER_CATCH()
//       }

//   // T can be statically converted to wstring_view
//   template<class T, typename std::enable_if<is_convertible_to_wstring_view<const T &>::value, int>::type = 0>
//   void log(source_loc loc, level::level_enum lvl, const T &msg)
//   {
//     bool log_enabled = should_log(lvl);
//     bool traceback_enabled = tracer_.enabled();
//     if (!log_enabled && !traceback_enabled)
//       {
// 	return;
//       }

//     SPDLOG_TRY
//       {
// 	memory_buf_t buf;
// 	details::os::wstr_to_utf8buf(msg, buf);
// 	details::log_msg log_msg(loc, name_, lvl, string_view_t(buf.data(), buf.size()));
// 	log_it_(log_msg, log_enabled, traceback_enabled, mpiRank_);
//       }
//     SPDLOG_LOGGER_CATCH()
//       }
// #endif // _WIN32
// #endif // SPDLOG_WCHAR_TO_UTF8_SUPPORT

  // return true logging is enabled for the given level.
  bool should_log(level::level_enum msg_level) const
  {
    return msg_level >= level_.load(std::memory_order_relaxed);
  }

  // return true if backtrace logging is enabled.
  bool should_backtrace() const
  {
    return tracer_.enabled();
  }

  void set_level(level::level_enum log_level)
  {
    level_.store(log_level);
  }

  level::level_enum level() const
  {
    return static_cast<level::level_enum>(level_.load(std::memory_order_relaxed));
  }

  const std::string &name() const
  {
    return name_;
  }


  // set formatting for the sinks in this logger.
  // each sink will get a separate instance of the formatter object.
  void set_formatter(std::unique_ptr<formatter> f)
  {
    for (auto it = sinks_.begin(); it != sinks_.end(); ++it)
      {
	if (std::next(it) == sinks_.end())
	  {
	    // last element - we can be move it.
	    (*it)->set_formatter(std::move(f));
	    break; // to prevent clang-tidy warning
	  }
	else
	  {
	    (*it)->set_formatter(f->clone());
	  }
      }
  }

  void set_pattern(std::string pattern, pattern_time_type time_type = pattern_time_type::local)
  {
    auto new_formatter = details::make_unique<pattern_formatter>(std::move(pattern), time_type);
    set_formatter(std::move(new_formatter));
  }

  // backtrace support.
  // efficiently store all debug/trace messages in a circular buffer until needed for debugging.
  void enable_backtrace(size_t n_messages)
  {
    tracer_.enable(n_messages);
  }

  void disable_backtrace()
  {
    tracer_.disable();
  }

  void dump_backtrace()
  {
    dump_backtrace_();
  }

  // flush functions
  void flush()
  {
    flush_();
  }

  void flush_on(level::level_enum log_level)
  {
    flush_level_.store(log_level);
  }

  level::level_enum flush_level() const
  {
    return static_cast<level::level_enum>(flush_level_.load(std::memory_order_relaxed));
  }

  // sinks
  const std::vector<sink_ptr> &sinks() const
  {
    return sinks_;
  }

  std::vector<sink_ptr> &sinks()
  {
    return sinks_;
  }

  // error handler
  void set_error_handler(err_handler handler)
  {
    custom_err_handler_ = std::move(handler);
  }

  // create new logger with same sinks and configuration.
  virtual std::shared_ptr<logger> clone(std::string logger_name)
  {
    auto cloned = std::make_shared<logger>(*this);
    cloned->name_ = std::move(logger_name);
    return cloned;
  }

protected:
  std::string name_;
  std::vector<sink_ptr> sinks_;
  spdlog::level_t level_{level::info};
  spdlog::level_t flush_level_{level::off};
  err_handler custom_err_handler_{nullptr};
  details::backtracer<bool> tracer_;

protected:
  // common implementation all log call after templated public api has been resolved
  template<typename FormatString, typename... Args>
  void log_(source_loc loc, level::level_enum lvl, const FormatString &fmt, Args&&...args)
  {
    bool log_enabled = should_log(lvl);
    bool traceback_enabled = tracer_.enabled();
    if (!log_enabled && !traceback_enabled)
    {
      return;
    }
    SPDLOG_TRY
    {
      memory_buf_t buf;
      fmt::format_to(buf, fmt, std::forward<Args>(args)...);
      details::log_msg log_msg(loc, name_, lvl, string_view_t(buf.data(), buf.size()));
      log_it_(log_msg, log_enabled, traceback_enabled);
    }
    SPDLOG_LOGGER_CATCH()
  }

  // log the given message (if the given log level is high enough),
  // and save backtrace (if backtrace is enabled).
  void log_it_(const details::log_msg &log_msg, bool log_enabled, bool traceback_enabled)
  {
    if (log_enabled)
    {
      sink_it_(log_msg);
    }
    if (traceback_enabled){
      tracer_.push_back(log_msg);
    }
  }

  virtual void sink_it_(const details::log_msg &msg)
  {
    for (auto &sink : sinks_)
    {
      if (sink->should_log(msg.level))
      {
	SPDLOG_TRY
	{
	  sink->log(msg);
	}
	SPDLOG_LOGGER_CATCH()
      }
    }

    if (should_flush_(msg)){
      flush_();
    }
  }

  virtual void flush_()
  {
    for (auto &sink : sinks_){
      SPDLOG_TRY
      {
	sink->flush();
      }
      SPDLOG_LOGGER_CATCH()
    }
  }

  void dump_backtrace_()
  {
    using details::log_msg;
    if (tracer_.enabled())
    {
      sink_it_(log_msg{name(), level::info, "****************** Backtrace Start ******************"});
      tracer_.foreach_pop([this](const log_msg &msg) { this->sink_it_(msg); });
      sink_it_(log_msg{name(), level::info, "****************** Backtrace End ********************"});
    }
  }

  bool should_flush_(const details::log_msg &msg)
  {
    auto flush_level = flush_level_.load(std::memory_order_relaxed);
    return (msg.level >= flush_level) && (msg.level != level::off);
  }

  // handle errors during logging.
  // default handler prints the error to stderr at max rate of 1 message/sec.
  void err_handler_(const std::string &msg)
  {
    if (custom_err_handler_)
      {
	custom_err_handler_(msg);
      }
    else
      {
	using std::chrono::system_clock;
	static std::mutex mutex;
	static std::chrono::system_clock::time_point last_report_time;
	static size_t err_counter = 0;
	std::lock_guard<std::mutex> lk{mutex};
	auto now = system_clock::now();
	err_counter++;
	if (now - last_report_time < std::chrono::seconds(1))
	  {
	    return;
	  }
	last_report_time = now;
	auto tm_time = details::os::localtime(system_clock::to_time_t(now));
	char date_buf[64];
	std::strftime(date_buf, sizeof(date_buf), "%Y-%m-%d %H:%M:%S", &tm_time);
#if defined(USING_R) && defined(R_R_H) // if in R environment
	REprintf("[*** LOG ERROR #%04zu ***] [%s] [%s] {%s}\n", err_counter, date_buf, name().c_str(), msg.c_str());
#else
	std::fprintf(stderr, "[*** LOG ERROR #%04zu ***] [%s] [%s] {%s}\n", err_counter, date_buf, name().c_str(), msg.c_str());
#endif
      }
  }
};

inline void swap(logger &a, logger &b)
{
  a.swap(b);
}

} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "logger-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_LOGGER_HPP_
