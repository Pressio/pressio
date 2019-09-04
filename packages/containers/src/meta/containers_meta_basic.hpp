/*
//@HEADER
// ************************************************************************
//
// containers_meta_basic.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_META_META_BASIC_HPP_
#define CONTAINERS_META_META_BASIC_HPP_

#include "containers_meta_has_communicator_typedef.hpp"
#include "containers_meta_has_data_map_typedef.hpp"
#include "containers_meta_has_global_ordinal_typedef.hpp"
#include "containers_meta_has_local_ordinal_typedef.hpp"
#include "containers_meta_has_ordinal_typedef.hpp"
#include "containers_meta_has_scalar_typedef.hpp"
#include "containers_meta_has_size_method.hpp"
#include "containers_meta_is_teuchos_rcp.hpp"

namespace pressio{ namespace containers{ namespace meta {

// template<typename T>
// struct remove_const: std::remove_const<T>{};

// template<typename T>
// struct remove_reference : std::remove_reference<T>{};

// template<typename T>
// struct remove_pointer : std::remove_pointer<T>{};

template<typename T>
struct is_arithmetic : std::is_arithmetic<T>{};

template<typename T>
struct is_integral: std::is_integral<T>{};

}}} // namespace pressio::containers::meta
#endif
