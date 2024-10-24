// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef PRESSIO_UTILS_LOGGER_SPDLOG_FORMATTER_HPP_
#define PRESSIO_UTILS_LOGGER_SPDLOG_FORMATTER_HPP_

// #include "./fmt/fmt.hpp"
//#include "./details/log_msg.hpp"

namespace spdlog {

class formatter
{
public:
    virtual ~formatter() = default;
    virtual void format(const details::log_msg &msg, memory_buf_t &dest) = 0;
    virtual std::unique_ptr<formatter> clone() const = 0;
};
} // namespace spdlog
#endif  // PRESSIO_UTILS_LOGGER_SPDLOG_FORMATTER_HPP_
