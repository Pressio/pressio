/*
//@HEADER
// ************************************************************************
//
//                          tuple_as_sequence.hpp
//                         dharma_new
//              Copyright (C) 2017 NTESS, LLC
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/

#ifndef SRC_TINYMPL_TUPLE_AS_SEQUENCE_HPP_
#define SRC_TINYMPL_TUPLE_AS_SEQUENCE_HPP_

#include <tuple>

#include "as_sequence.hpp"
#include "sequence.hpp"


namespace tinympl {

template <typename... Args>
struct as_sequence<std::tuple<Args...>>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = std::tuple<Ts...>;
};

template <typename... Args>
struct as_sequence<const std::tuple<Args...>>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const std::tuple<Ts...>;
};

template <typename... Args>
struct as_sequence<const volatile std::tuple<Args...>>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const volatile std::tuple<Ts...>;
};

template <typename... Args>
struct as_sequence<volatile std::tuple<Args...>>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = volatile std::tuple<Ts...>;
};

template <typename... Args>
struct as_sequence<std::tuple<Args...>&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = std::tuple<Ts...>&;
};

template <typename... Args>
struct as_sequence<const std::tuple<Args...>&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const std::tuple<Ts...>&;
};

template <typename... Args>
struct as_sequence<const volatile std::tuple<Args...>&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const volatile std::tuple<Ts...>&;
};

template <typename... Args>
struct as_sequence<volatile std::tuple<Args...>&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = volatile std::tuple<Ts...>&;
};

template <typename... Args>
struct as_sequence<std::tuple<Args...>&&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = std::tuple<Ts...>&&;
};

template <typename... Args>
struct as_sequence<const std::tuple<Args...>&&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const std::tuple<Ts...>&&;
};

template <typename... Args>
struct as_sequence<const volatile std::tuple<Args...>&&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = const volatile std::tuple<Ts...>&&;
};

template <typename... Args>
struct as_sequence<volatile std::tuple<Args...>&&>
{
  typedef sequence<Args...> type;
  template <class... Ts> using rebind = volatile std::tuple<Ts...>&&;
};


} // end namespace tinympl


#endif /* SRC_TINYMPL_TUPLE_AS_SEQUENCE_HPP_ */
