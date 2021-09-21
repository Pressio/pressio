// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_ASYNC_LOGGER_HPP_
#define UTILS_LOGGER_SPDLOG_ASYNC_LOGGER_HPP_

// Fast asynchronous logger.
// Uses pre allocated queue.
// Creates a single back thread to pop messages from the queue and log them.
//
// Upon each log write the logger:
//    1. Checks if its log level is enough to log the message
//    2. Push a new copy of the message to a queue (or block the caller until
//    space is available in the queue)
// Upon destruction, logs all remaining messages in the queue before
// destructing..

#include "./logger.hpp"
#include "./sinks/sink.hpp"
#include "./details/thread_pool.hpp"

#include <memory>
#include <string>

namespace spdlog {

// Async overflow policy - block by default.
enum class async_overflow_policy
  {
    block,         // Block until message can be enqueued
    overrun_oldest // Discard oldest message in the queue if full when trying to
                   // add new item.
  };

namespace details {
class thread_pool;
}

class async_logger final : public std::enable_shared_from_this<async_logger>, public logger
{
  friend class details::thread_pool;

public:
  template<typename It>
  async_logger(std::string logger_name,
	       It begin, It end,
	       std::weak_ptr<details::thread_pool> tp,
	       async_overflow_policy overflow_policy = async_overflow_policy::block)
    : logger(std::move(logger_name), begin, end)
    , thread_pool_(std::move(tp))
    , overflow_policy_(overflow_policy)
  {}

  async_logger(std::string logger_name,
	       sinks_init_list sinks_list,
	       std::weak_ptr<details::thread_pool> tp,
	       async_overflow_policy overflow_policy = async_overflow_policy::block)
    : async_logger(std::move(logger_name),
		   sinks_list.begin(),
		   sinks_list.end(),
		   std::move(tp),
		   overflow_policy)
  {}

  async_logger(std::string logger_name,
	       sink_ptr single_sink,
	       std::weak_ptr<details::thread_pool> tp,
	       async_overflow_policy overflow_policy = async_overflow_policy::block)
    : async_logger(std::move(logger_name), {std::move(single_sink)}, std::move(tp), overflow_policy)
  {}

  std::shared_ptr<logger> clone(std::string new_name) override
  {
    auto cloned = std::make_shared<spdlog::async_logger>(*this);
    cloned->name_ = std::move(new_name);
    return cloned;
  }

protected:
  void sink_it_(const details::log_msg &msg) override
  {
    if (auto pool_ptr = thread_pool_.lock())
      {
        pool_ptr->post_log(shared_from_this(), msg, overflow_policy_);
      }
    else
      {
        throw_spdlog_ex("async log: thread pool doesn't exist anymore");
      }
  }

  void flush_() override
  {
    if (auto pool_ptr = thread_pool_.lock())
      {
        pool_ptr->post_flush(shared_from_this(), overflow_policy_);
      }
    else
      {
        throw_spdlog_ex("async flush: thread pool doesn't exist anymore");
      }
  }

  void backend_sink_it_(const details::log_msg &incoming_log_msg)
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

    if (should_flush_(msg))
      {
        backend_flush_();
      }
  }

  void backend_flush_()
  {
    for (auto &sink : sinks_)
      {
        SPDLOG_TRY
	  {
            sink->flush();
	  }
        SPDLOG_LOGGER_CATCH()
	  }
  }

private:
  std::weak_ptr<details::thread_pool> thread_pool_;
  async_overflow_policy overflow_policy_;
};
} // namespace spdlog

// #ifdef SPDLOG_HEADER_ONLY
// #include "async_logger-inl.hpp"
// #endif
#endif  // UTILS_LOGGER_SPDLOG_ASYNC_LOGGER_HPP_
