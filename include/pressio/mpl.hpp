/*
//@HEADER
// ************************************************************************
//
// mpl.hpp
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

#ifndef PRESSIO_MPL_HPP_
#define PRESSIO_MPL_HPP_

#include "macros.hpp"

#include "./mpl/mpl_ConfigDefs.hpp"

// some will change/disappear once we move to C++14
#include "./mpl/identity.hpp"
#include "./mpl/enable_if_t.hpp"
#include "./mpl/conditional_t.hpp"
#include "./mpl/void_t.hpp"
#include "./mpl/not_same.hpp"
#include "./mpl/remove_cvref.hpp"
#include "./mpl/remove_reference.hpp"
#include "./mpl/not_void.hpp"
#include "./mpl/all_of.hpp"
#include "./mpl/none_of.hpp"
#include "./mpl/any_of.hpp"
#include "./mpl/at.hpp"
#include "./mpl/size.hpp"

#include "./mpl/detection_idiom.hpp"
#include "./mpl/is_subscriptable_as.hpp"
#include "./mpl/is_default_constructible.hpp"
#include "./mpl/is_std_complex.hpp"
#include "./mpl/is_std_shared_ptr.hpp"
#include "./mpl/is_std_unique_ptr.hpp"
#include "./mpl/is_std_shared_ptr.hpp"
#include "./mpl/publicly_inherits_from.hpp"

// #include "./mpl/mpl_non_variadic.hpp"
#include "./mpl/mpl_variadic.hpp"

#endif
