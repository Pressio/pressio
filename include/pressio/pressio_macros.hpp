/*
//@HEADER
// ************************************************************************
//
// macros.hpp
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

#ifndef PRESSIO_MACROS_HPP_
#define PRESSIO_MACROS_HPP_

#include "pressio/ops_macros.hpp"

// ----------------------------------------
// logging macros
// ----------------------------------------
#define PRESSIO_LOG_LEVEL_TRACE	    0
#define PRESSIO_LOG_LEVEL_DEBUG	    1
#define PRESSIO_LOG_LEVEL_INFO	    2
#define PRESSIO_LOG_LEVEL_WARN	    3
#define PRESSIO_LOG_LEVEL_ERROR	    4
#define PRESSIO_LOG_LEVEL_CRITICAL  5
#define PRESSIO_LOG_LEVEL_OFF	    6

// if we are in debug mode, enable debug prints by default
#if !defined NDEBUG && !defined PRESSIO_ENABLE_DEBUG_PRINT
#define PRESSIO_ENABLE_DEBUG_PRINT
#endif

#if defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
// if DEBUG_PRINT is on but MIN_LEVEL off, then set min level to trace
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_TRACE
#elif !defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
// if DEBUG_PRINT=off and MIN_LEVEL=off, then set logging off
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_OFF
#else
//DEBUG_PRINT is off and LOG_ACTIVE_MIN_LEVEL=on, nothing to do
#endif

#endif
