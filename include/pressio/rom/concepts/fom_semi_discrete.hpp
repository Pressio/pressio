/*
//@HEADER
// ************************************************************************
//
// rom_fom_system_continuous_time.hpp
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

#ifndef ROM_CONCEPTS_FOM_SEMI_DISCRETE_HPP_
#define ROM_CONCEPTS_FOM_SEMI_DISCRETE_HPP_

#include "helpers.hpp"

#ifdef PRESSIO_ENABLE_CXX20

// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace rom{

template <class T>
concept SemiDiscreteFom =
  /* must have nested aliases */
  requires(){
    typename T::time_type;
    typename T::state_type;
    typename T::right_hand_side_type;
  }
  /*
    requirements on the nested aliases
  */
  && std::regular<typename T::time_type>
  && std::totally_ordered<typename T::time_type>
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::right_hand_side_type>
  && all_have_traits_and_same_scalar<
    typename T::state_type, typename T::right_hand_side_type>::value
  && Traits<typename T::state_type>::rank == 1
  && Traits<typename T::right_hand_side_type>::rank == 1
  && std::convertible_to<
    typename T::time_type, scalar_trait_t<typename T::state_type>>
  /*
    compound and nested requirements last
  */
  && requires(const T & A,
	      const typename T::state_type     & state,
	      const typename T::time_type      & evalTime,
	      typename T::right_hand_side_type & rhs)
  {
    { A.createRightHandSide() } -> std::same_as<typename T::right_hand_side_type>;
    { A.rightHandSide(state, evalTime, rhs) } -> std::same_as<void>;
  };

}} // end namespace pressio::rom







/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else

namespace pressio{ namespace rom{

template<class T, class enable = void>
struct SemiDiscreteFom : std::false_type{};

template<class T>
struct SemiDiscreteFom<
  T,
  mpl::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_right_hand_side_typedef<T>::value
    //
    // regular and totally ordered are in C++20
    // && std::regular<typename T::time_type>
    // && std::totally_ordered<typename T::time_type>
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::right_hand_side_type>::value
    && all_have_traits_and_same_scalar<
       typename T::state_type, typename T::right_hand_side_type>::value
    && Traits<typename T::state_type>::rank == 1
    && Traits<typename T::right_hand_side_type>::rank == 1
    && std::is_convertible<
	 typename T::time_type, scalar_trait_t<typename T::state_type>
       >::value
    //
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::right_hand_side_type >::value
    && ::pressio::rom::has_const_rhs_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::right_hand_side_type>::value
   >
  > : std::true_type{};

}} // end namespace pressio::rom

#endif

#endif  // ROM_CONCEPTS_FOM_SEMI_DISCRETE_HPP_
