/*
//@HEADER
// ************************************************************************
//
// pressio_utils.hpp
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

#ifndef PRESSIO_UTILS_HPP_
#define PRESSIO_UTILS_HPP_

#include "pressio_mpl.hpp"

#include "utils/src/utils_crtp_helper.hpp"
#include "utils/src/utils_static_constants.hpp"
#include "utils/src/utils_empty.hpp"
#include "utils/src/utils_instance_or_reference_wrapper.hpp"
#include "utils/src/utils_read_ascii_matrix_std_vec_vec.hpp"
#include "utils/src/utils_set_stream_precision.hpp"

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
#include "utils/src/utils_teuchos_performance_monitor.hpp"
#endif

#include "utils/src/io/utils_colorize_print.hpp"
#include "utils/src/io/utils_print_helper.hpp"

/* headers needed for logging */

// fmt is needed
#include "utils/src/logger/fmt/fmt.hpp"

#define PRESSIO_LOG_LEVEL_TRACE		0
#define PRESSIO_LOG_LEVEL_DEBUG		1
#define PRESSIO_LOG_LEVEL_INFO		2
#define PRESSIO_LOG_LEVEL_WARN		3
#define PRESSIO_LOG_LEVEL_ERROR		4
#define PRESSIO_LOG_LEVEL_CRITICAL	5
#define PRESSIO_LOG_LEVEL_OFF		6

#if defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_TRACE
#elif !defined(PRESSIO_ENABLE_DEBUG_PRINT) && !defined(PRESSIO_LOG_ACTIVE_MIN_LEVEL)
#define PRESSIO_LOG_ACTIVE_MIN_LEVEL	PRESSIO_LOG_LEVEL_OFF
#endif

#include "utils/src/logger/utils_logger_enums.hpp"
#include "utils/src/logger/utils_logger.hpp"

#endif
