// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_SINKS_UNUSED_BASIC_FILE_SINK_INL_HPP_
#define UTILS_LOGGER_SPDLOG_SINKS_UNUSED_BASIC_FILE_SINK_INL_HPP_

// #ifndef SPDLOG_HEADER_ONLY
// #include "../sinks/basic_file_sink.hpp"
// #endif

// #include "../common.hpp"
// #include "../details/os.hpp"

namespace spdlog {
namespace sinks {

// template<typename Mutex>
// SPDLOG_INLINE basic_file_sink<Mutex>::basic_file_sink(const filename_t &filename, bool truncate)
// {
//     file_helper_.open(filename, truncate);
// }

// template<typename Mutex>
// SPDLOG_INLINE const filename_t &basic_file_sink<Mutex>::filename() const
// {
//     return file_helper_.filename();
// }

// template<typename Mutex>
// SPDLOG_INLINE void basic_file_sink<Mutex>::sink_it_(const details::log_msg &msg)
// {
//     memory_buf_t formatted;
//     base_sink<Mutex>::formatter_->format(msg, formatted);
//     file_helper_.write(formatted);
// }

// template<typename Mutex>
// SPDLOG_INLINE void basic_file_sink<Mutex>::flush_()
// {
//     file_helper_.flush();
// }

} // namespace sinks
} // namespace spdlog
#endif  // UTILS_LOGGER_SPDLOG_SINKS_UNUSED_BASIC_FILE_SINK_INL_HPP_
