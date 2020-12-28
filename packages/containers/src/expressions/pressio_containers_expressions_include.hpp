/*
//@HEADER
// ************************************************************************
//
// pressio_containers_expressions_include.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_PRESSIO_CONTAINERS_EXPRESSIONS_INCLUDE_HPP_
#define CONTAINERS_EXPRESSIONS_PRESSIO_CONTAINERS_EXPRESSIONS_INCLUDE_HPP_

#include "./span/containers_span_classes.hpp"
#include "./span/containers_span_traits.hpp"
#include "./span/containers_span_function.hpp"
#include "./span/containers_is_expression_span.hpp"

#include "./subspan/containers_subspan_classes.hpp"
#include "./subspan/containers_subspan_traits.hpp"
#include "./subspan/containers_subspan_function.hpp"
#include "./subspan/containers_is_expression_subspan.hpp"

#include "./diag/containers_diag_classes.hpp"
#include "./diag/containers_diag_traits.hpp"
#include "./diag/containers_diag_function.hpp"
#include "./diag/containers_is_expression_diag.hpp"

#include "./asdiagonalmatrix/containers_asdiagonalmatrix_classes.hpp"
#include "./asdiagonalmatrix/containers_asdiagonalmatrix_traits.hpp"
#include "./asdiagonalmatrix/containers_asdiagonalmatrix_function.hpp"
#include "./asdiagonalmatrix/containers_is_expression_asdiagonalmatrix.hpp"

#include "containers_is_expression.hpp"

#endif  // CONTAINERS_EXPRESSIONS_PRESSIO_CONTAINERS_EXPRESSIONS_INCLUDE_HPP_
