// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_FORMATTER_H_
#define UTILS_LOGGER_SPDLOG_FORMATTER_H_

#include "./fmt/fmt.h"
#include "./details/log_msg.h"

namespace spdlog {

class formatter
{
public:
    virtual ~formatter() = default;
    virtual void format(const details::log_msg &msg, memory_buf_t &dest) = 0;
    virtual std::unique_ptr<formatter> clone() const = 0;
};
} // namespace spdlog
#endif  // UTILS_LOGGER_SPDLOG_FORMATTER_H_
