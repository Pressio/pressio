
#ifndef UTILS_UTILS_LOGGER_ENUMS_HPP_
#define UTILS_UTILS_LOGGER_ENUMS_HPP_

namespace pressio{

enum class logto{terminal, file, fileAndTerminal, terminalAndFile};

namespace log{
enum class level
{
  trace		= PRESSIO_LOG_LEVEL_TRACE,
  debug		= PRESSIO_LOG_LEVEL_DEBUG,
  info		= PRESSIO_LOG_LEVEL_INFO,
  warn		= PRESSIO_LOG_LEVEL_WARN,
  err		= PRESSIO_LOG_LEVEL_ERROR,
  critical	= PRESSIO_LOG_LEVEL_CRITICAL,
  off		= PRESSIO_LOG_LEVEL_OFF,
  n_levels
};
}// end namespace pressio::log
}// end namespace pressio

#endif  // UTILS_UTILS_LOGGER_ENUMS_HPP_
