// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef UTILS_LOGGER_SPDLOG_PATTERN_FORMATTER_HPP_
#define UTILS_LOGGER_SPDLOG_PATTERN_FORMATTER_HPP_

// #include "./common.hpp"
// #include "./details/log_msg.hpp"
//#include "./details/os.hpp"
//#include "./formatter.hpp"
//#include "./details/fmt_helper.hpp"
//#include "./details/log_msg.hpp"
//#include "./details/os.hpp"
// #include "./fmt/fmt.hpp"
//#include "./formatter.hpp"

// #include <chrono>
// #include <ctime>
// #include <memory>
// #include <string>
// #include <vector>
#include <unordered_map>

namespace spdlog {
namespace details {

// padding information.
struct padding_info
{
  enum class pad_side
    {
      left,
      right,
      center
    };

  padding_info() = default;
  padding_info(size_t width, padding_info::pad_side side, bool truncate)
    : width_(width)
    , side_(side)
    , truncate_(truncate)
    , enabled_(true)
  {}

  bool enabled() const
  {
    return enabled_;
  }
  size_t width_ = 0;
  pad_side side_ = pad_side::left;
  bool truncate_ = false;
  bool enabled_ = false;
};

class flag_formatter
{
public:
  explicit flag_formatter(padding_info padinfo)
    : padinfo_(padinfo)
  {}
  flag_formatter() = default;
  virtual ~flag_formatter() = default;
  virtual void format(const details::log_msg &msg, const std::tm &tm_time, memory_buf_t &dest) = 0;

protected:
  padding_info padinfo_;
};

} // namespace details

class custom_flag_formatter : public details::flag_formatter
{
public:
  virtual std::unique_ptr<custom_flag_formatter> clone() const = 0;

  void set_padding_info(details::padding_info padding)
  {
    flag_formatter::padinfo_ = padding;
  }
};

namespace details{
class scoped_padder;
class full_formatter;
class aggregate_formatter;
struct null_scoped_padder;
}

class pattern_formatter final : public formatter
{

public:
  using custom_flags = std::unordered_map<char, std::unique_ptr<custom_flag_formatter>>;

private:
  std::string pattern_;
  std::string eol_;
  pattern_time_type pattern_time_type_;
  std::tm cached_tm_;
  std::chrono::seconds last_log_secs_;
  std::vector<std::unique_ptr<details::flag_formatter>> formatters_;
  custom_flags custom_handlers_;

public:
  explicit pattern_formatter(std::string pattern,
			     pattern_time_type time_type = pattern_time_type::local,
			     std::string eol = spdlog::details::os::default_eol,
			     custom_flags custom_user_flags = custom_flags());

  // use default pattern is not given
  explicit pattern_formatter(pattern_time_type time_type = pattern_time_type::local,
			     std::string eol = spdlog::details::os::default_eol);

  pattern_formatter(const pattern_formatter &other) = delete;
  pattern_formatter &operator=(const pattern_formatter &other) = delete;

public:
  std::unique_ptr<formatter> clone() const override
  {
    custom_flags cloned_custom_formatters;
    for (auto &it : custom_handlers_)
      {
        cloned_custom_formatters[it.first] = it.second->clone();
      }
    return details::make_unique<pattern_formatter>(pattern_, pattern_time_type_, eol_, std::move(cloned_custom_formatters));
  }

  void format(const details::log_msg &msg, memory_buf_t &dest) override
  {
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(msg.time.time_since_epoch());
    if (secs != last_log_secs_)
      {
        cached_tm_ = get_time_(msg);
        last_log_secs_ = secs;
      }

    for (auto &f : formatters_)
      {
        f->format(msg, cached_tm_, dest);
      }
    details::fmt_helper::append_string_view(eol_, dest);
  }

  template<typename T, typename... Args>
  pattern_formatter &add_flag(char flag, Args&&...args)
  {
    custom_handlers_[flag] = details::make_unique<T>(std::forward<Args>(args)...);
    return *this;
  }

  void set_pattern(std::string pattern)
  {
    pattern_ = std::move(pattern);
    compile_pattern_(pattern_);
  }

private:
  std::tm get_time_(const details::log_msg &msg)
  {
    if (pattern_time_type_ == pattern_time_type::local)
      {
	return details::os::localtime(log_clock::to_time_t(msg.time));
      }
    return details::os::gmtime(log_clock::to_time_t(msg.time));
  }

  // this one is needed if we want to change the pattern
  // from the default one
  template<typename Padder>
  void handle_flag_(char flag, details::padding_info padding);

  // Extract given pad spec (e.g. %8X)
  // Advance the given it pass the end of the padding spec found (if any)
  // Return padding.
  static details::padding_info handle_padspec_(std::string::const_iterator &it, std::string::const_iterator end)
  {
    using details::padding_info;
    using details::scoped_padder;
    const size_t max_width = 64;
    if (it == end)
      {
	return padding_info{};
      }

    padding_info::pad_side side;
    switch (*it)
      {
      case '-':
	side = padding_info::pad_side::right;
	++it;
	break;
      case '=':
	side = padding_info::pad_side::center;
	++it;
	break;
      default:
	side = details::padding_info::pad_side::left;
	break;
      }

    if (it == end || !std::isdigit(static_cast<unsigned char>(*it)))
      {
	return padding_info{}; // no padding if no digit found here
      }

    auto width = static_cast<size_t>(*it) - '0';
    for (++it; it != end && std::isdigit(static_cast<unsigned char>(*it)); ++it)
      {
	auto digit = static_cast<size_t>(*it) - '0';
	width = width * 10 + digit;
      }

    // search for the optional truncate marker '!'
    bool truncate;
    if (it != end && *it == '!')
      {
	truncate = true;
	++it;
      }
    else
      {
	truncate = false;
      }

    return details::padding_info{std::min<size_t>(width, max_width), side, truncate};
  }

  void compile_pattern_(const std::string &pattern);

};
} // namespace spdlog

#ifdef SPDLOG_HEADER_ONLY
#include "pattern_formatter-inl.hpp"
#endif

#endif  // UTILS_LOGGER_SPDLOG_PATTERN_FORMATTER_HPP_
