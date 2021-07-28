// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_SINKS_BASIC_FILE_SINK_HPP_
#define UTILS_LOGGER_SPDLOG_SINKS_BASIC_FILE_SINK_HPP_

#include "../details/file_helper.hpp"
//#include "../details/null_mutex.hpp"
#include "../sinks/base_sink.hpp"
//#include "../details/synchronous_factory.hpp"
// #include <mutex>
// #include <string>

namespace spdlog { namespace sinks {

/* Trivial file sink with single file as target */
template<typename Mutex>
class basic_file_sink final : public base_sink<Mutex>
{
public:
  explicit basic_file_sink(const filename_t &filename, bool truncate = false)
  {
    file_helper_.open(filename, truncate);
  }

  const filename_t &filename() const
  {
    return file_helper_.filename();
  }


protected:
  void sink_it_(const details::log_msg &msg) override
  {
    memory_buf_t formatted;
    base_sink<Mutex>::formatter_->format(msg, formatted);
    file_helper_.write(formatted);
  }

  void flush_() override
  {
    file_helper_.flush();
  }


private:
  details::file_helper file_helper_;
};

using basic_file_sink_mt = basic_file_sink<std::mutex>;
using basic_file_sink_st = basic_file_sink<details::null_mutex>;

} // namespace sinks

// //
// // factory functions
// //
// template<typename Factory = spdlog::synchronous_factory>
// inline std::shared_ptr<logger> basic_logger_mt(const std::string &logger_name, const filename_t &filename, bool truncate = false)
// {
//   return Factory::template create<sinks::basic_file_sink_mt>(logger_name, filename, truncate);
// }

// template<typename Factory = spdlog::synchronous_factory>
// inline std::shared_ptr<logger> basic_logger_st(const std::string &logger_name, const filename_t &filename, bool truncate = false)
// {
//   return Factory::template create<sinks::basic_file_sink_st>(logger_name, filename, truncate);
// }

} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "basic_file_sink-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_SINKS_BASIC_FILE_SINK_HPP_
