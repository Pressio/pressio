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

#ifndef ROM_CONSTRAINTS_ROM_SEMI_DISCRETE_FOM_CONCEPT_HPP_
#define ROM_CONSTRAINTS_ROM_SEMI_DISCRETE_FOM_CONCEPT_HPP_

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20

template <class T>
concept SemiDiscreteFom =
      std::regular<typename T::time_type>
   && std::totally_ordered<typename T::time_type>
   && std::copy_constructible<typename T::state_type>
   && std::copy_constructible<typename T::right_hand_side_type>
   && std::same_as<
       typename pressio::Traits<typename T::state_type>::scalar_type,
       typename pressio::Traits<typename T::right_hand_side_type>::scalar_type>
   && std::convertible_to<
       typename T::time_type,
       typename pressio::Traits<typename T::state_type>::scalar_type>
   && requires(const T & A,
               const typename T::state_type & s,
               const typename T::time_type & t,
               typename T::right_hand_side_type & r)
   {
    { A.createRightHandSide()  } -> std::same_as<typename T::right_hand_side_type>;
    { A.rightHandSide(s, t, r) } -> std::same_as<void>;
  };

#else

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
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::right_hand_side_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::right_hand_side_type
      >::value
    && std::is_convertible<
      typename T::time_type,
      typename ::pressio::Traits<typename T::state_type>::scalar_type
      >::value
    //
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::right_hand_side_type >::value
    && ::pressio::ode::has_const_rhs_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::right_hand_side_type>::value
   >
  > : std::true_type{};

#endif

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
