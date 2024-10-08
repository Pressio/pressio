
// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_DETAILS_BACKTRACER_HPP_
#define UTILS_LOGGER_SPDLOG_DETAILS_BACKTRACER_HPP_

//#include "log_msg_buffer.hpp"
//#include "circular_q.hpp"
//#include <atomic>
#include <mutex>
#include <functional>

// Store log messages in circular buffer.
// Useful for storing debug data in case of error/warning happens.

namespace spdlog {
namespace details {

template<typename T = bool>
class backtracer
{
  mutable std::mutex mutex_;
  std::atomic<T> enabled_{false};
  circular_q<log_msg_buffer> messages_;

public:
  backtracer() = default;

  backtracer(const backtracer &other)
  {
    std::lock_guard<std::mutex> lock(other.mutex_);
    enabled_ = other.enabled();
    messages_ = other.messages_;
  }

  backtracer(backtracer &&other) SPDLOG_NOEXCEPT
  {
    std::lock_guard<std::mutex> lock(other.mutex_);
    enabled_ = other.enabled();
    messages_ = std::move(other.messages_);
  }

  backtracer &operator=(backtracer other)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    enabled_ = other.enabled();
    messages_ = std::move(other.messages_);
    return *this;
  }

  void enable(size_t size)
  {
    std::lock_guard<std::mutex> lock{mutex_};
    enabled_.store(true, std::memory_order_relaxed);
    messages_ = circular_q<log_msg_buffer>{size};
  }

  void disable()
  {
    std::lock_guard<std::mutex> lock{mutex_};
    enabled_.store(false, std::memory_order_relaxed);
  }

  bool enabled() const
  {
    return enabled_.load(std::memory_order_relaxed);
  }

  void push_back(const log_msg &msg)
  {
    std::lock_guard<std::mutex> lock{mutex_};
    messages_.push_back(log_msg_buffer{msg});
  }

  // pop all items in the q and apply the given fun on each of them.
  void foreach_pop(std::function<void(const details::log_msg &)> fun)
  {
    std::lock_guard<std::mutex> lock{mutex_};
    while (!messages_.empty())
      {
        auto &front_msg = messages_.front();
        fun(front_msg);
        messages_.pop_front();
      }
  }
};

} // namespace details
} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "backtracer-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_DETAILS_BACKTRACER_HPP_
