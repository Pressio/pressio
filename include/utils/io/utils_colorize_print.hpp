/*
//@HEADER
// ************************************************************************
//
// utils_colorize_print.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef UTILS_IO_UTILS_COLORIZE_PRINT_HPP_
#define UTILS_IO_UTILS_COLORIZE_PRINT_HPP_

#include <iostream>
#include <unistd.h>

namespace pressio{ namespace utils{ namespace io{

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

}}} // namespace pressio::utils::io


#endif  // UTILS_IO_UTILS_COLORIZE_PRINT_HPP_
