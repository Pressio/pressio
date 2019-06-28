
#ifndef ALGEBRA_COLORIZE_HPP_
#define ALGEBRA_COLORIZE_HPP_

#include <unistd.h>
#include <iostream>
#include <cstdio>

namespace rompp{ namespace algebra{ namespace io{

namespace impl{

inline
FILE* get_standard_stream(const std::ostream& stream){
  if (&stream == &std::cout)
    return stdout;
  else if ((&stream == &std::cerr) || (&stream == &std::clog))
    return stderr;
  return 0;
}

//! test if a `std::ostream` object refers to a terminal.
inline
bool is_atty(const std::ostream& stream){
  FILE* std_stream = get_standard_stream(stream);
  // assume it's not a tty if standard stream not detected
  if (!std_stream) return false;
  return ::isatty(fileno(std_stream));
}

inline bool is_colorized(std::ostream& stream){
  return is_atty(stream);
}

}//end namepsace utils::io::impl



// reset is needed after any command to reset default color
inline
std::string reset(){ return impl::is_colorized(std::cout) ? "\033[00m" : ""; }

// features
inline
std::string bold(){ return impl::is_colorized(std::cout) ? "\033[1m" : ""; }

inline
std::string dark(){ return impl::is_colorized(std::cout) ? "\033[2m" : ""; }

inline
std::string underline(){ return impl::is_colorized(std::cout) ? "\033[4m" : ""; }

inline
std::string blink(){ return impl::is_colorized(std::cout) ? "\033[5m" : ""; }


// colors
inline
std::string grey(){ return impl::is_colorized(std::cout) ? "\033[30m" : ""; }

inline
std::string red(){ return impl::is_colorized(std::cout) ? "\033[31m" : ""; }

inline
std::string green(){ return impl::is_colorized(std::cout) ? "\033[32m" : ""; }

inline
std::string yellow(){ return impl::is_colorized(std::cout) ? "\033[33m" : ""; }

inline
std::string blue(){ return impl::is_colorized(std::cout) ? "\033[34m" : ""; }

inline
std::string magenta(){ return impl::is_colorized(std::cout) ? "\033[35m" : ""; }

inline
std::string cyan(){ return impl::is_colorized(std::cout) ? "\033[36m" : ""; }

inline
std::string white(){ return impl::is_colorized(std::cout) ? "\033[37m" : ""; }


// background colors
inline
std::string bg_grey(){ return impl::is_colorized(std::cout) ? "\033[40m" : ""; }

inline
std::string bg_red(){ return impl::is_colorized(std::cout) ? "\033[41m" : ""; }

inline
std::string bg_green(){ return impl::is_colorized(std::cout) ? "\033[42m" : ""; }

inline
std::string bg_yellow(){ return impl::is_colorized(std::cout) ? "\033[43m" : ""; }

inline
std::string bg_blue(){ return impl::is_colorized(std::cout) ? "\033[44m" : ""; }

inline
std::string bg_magenta(){ return impl::is_colorized(std::cout) ? "\033[45m" : ""; }

inline
std::string bg_cyan(){ return impl::is_colorized(std::cout) ? "\033[46m" : ""; }

inline
std::string bg_white(){ return impl::is_colorized(std::cout) ? "\033[47m" : ""; }

}}} // namespace rompp::utils::io

#endif








// inline
// std::ostream& reset(std::ostream& stream){
//   if (impl::is_colorized(stream))stream << "\033[00m";
//   return stream;
// }

// inline
// std::ostream& bold(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[1m";
//   return stream;
// }

// inline
// std::ostream& dark(std::ostream& stream)
// {
//   if (impl::is_colorized(stream)) stream << "\033[2m";
//   return stream;
// }

// inline
// std::ostream& underline(std::ostream& stream)
// {
//   if (impl::is_colorized(stream)) stream << "\033[4m";
//   return stream;
// }

// inline
// std::ostream& blink(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[5m";
//   return stream;
// }

// // inline
// // std::ostream& reverse(std::ostream& stream){
// //   if (impl::is_colorized(stream)) stream << "\033[7m";
// //   return stream;
// // }

// // inline
// // std::ostream& concealed(std::ostream& stream){
// //   if (impl::is_colorized(stream)) stream << "\033[8m";
// //   return stream;
// // }

// inline
// std::ostream& grey(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[30m";
//   return stream;
// }

// inline
// std::ostream& red(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[31m";
//   return stream;
// }

// inline
// std::ostream& green(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[32m";
//   return stream;
// }

// inline
// std::ostream& yellow(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[33m";
//   return stream;
// }

// inline
// std::ostream& blue(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[34m";
//   return stream;
// }

// inline
// std::ostream& magenta(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[35m";
//   return stream;
// }

// inline
// std::ostream& cyan(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[36m";
//   return stream;
// }

// inline
// std::ostream& white(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[37m";
//   return stream;
// }

// inline
// std::ostream& on_grey(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[40m";
//   return stream;
// }

// inline
// std::ostream& on_red(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[41m";
//   return stream;
// }

// inline
// std::ostream& on_green(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[42m";
//   return stream;
// }

// inline
// std::ostream& on_yellow(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[43m";
//   return stream;
// }

// inline
// std::ostream& on_blue(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[44m";
//   return stream;
// }

// inline
// std::ostream& on_magenta(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[45m";
//   return stream;
// }

// inline
// std::ostream& on_cyan(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[46m";
//   return stream;
// }

// inline
// std::ostream& on_white(std::ostream& stream){
//   if (impl::is_colorized(stream)) stream << "\033[47m";
//   return stream;
// }
