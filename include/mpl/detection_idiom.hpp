/*
//@HEADER
// ************************************************************************
//
// detection_idiom.hpp
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

#ifndef MPL_DETECTION_IDIOM_HPP_
#define MPL_DETECTION_IDIOM_HPP_

namespace pressio{ namespace mpl{
  
struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

// already defined in containers_meta_basic
// template <typename...>
// using my_void_t = void;

template <class Default,
	  class AlwaysVoid,
	  template <class...> class Op,
	  class...>
struct detector {
  constexpr static auto value = false;  
  using type = Default;
};

template <class Default,
	  template <class...> class Op,
	  class... Args>
struct detector<Default, 
                ::pressio::mpl::void_t<Op<Args...>>, Op, Args...> {
  constexpr static auto value = true;
  using type = Op<Args...>;
};

//================================================
  
template <template <class...> class Op, class... Args>
using is_detected = detector<nonesuch, void, Op, Args...>;

template <template <class...> class Op, class... Args>
using detected_t = typename detector<nonesuch, void, Op, Args...>::type;
  
template <class T, template<class...> class Op, class... Args>
using detected_or = detector<T, void, Op, Args...>;  

template <class T, template<class...> class Op, class... Args>
using detected_or_t = typename detected_or<T, Op, Args...>::type;

template <class T, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<T, detected_t<Op, Args...>>;
  

}}//end namespace pressio::mpl
#endif  // MPL_DETECTION_IDIOM_HPP_
