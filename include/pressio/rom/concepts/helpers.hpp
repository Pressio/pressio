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

#ifndef ROM_CONSTRAINTS_ROM_CONCEPTS_HELPERS_HPP_
#define ROM_CONSTRAINTS_ROM_CONCEPTS_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace concepts{ namespace impl{

template<class FomSystemType, class OperandType>
struct fully_discrete_fom_jacobian_action{
  using type =
    decltype
    (std::declval<FomSystemType const>().createResultOfDiscreteTimeJacobianActionOn
     (
      std::declval<OperandType const &>()
     )
    );
};

template<class ...Args>
using fully_discrete_fom_jacobian_action_t =
  typename fully_discrete_fom_jacobian_action<Args...>::type;

// -----------------------------------------------------------------------------

template<class FomSystemType, class TrialSubspaceType>
struct fom_jacobian_action_on_trial_space{
  using type =
    decltype
    (std::declval<FomSystemType const>().createApplyJacobianResult
     (
      std::declval<typename TrialSubspaceType::basis_matrix_type const &>()
     )
    );
};

template<class ...Args>
using fom_jacobian_action_on_trial_space_t =
  typename fom_jacobian_action_on_trial_space<Args...>::type;

// -----------------------------------------------------------------------------

template<class FomSystemType, class OperandType>
struct fom_jacobian_action{
  using type =
    decltype
    (std::declval<FomSystemType const>().createApplyJacobianResult
     (
      std::declval<OperandType const &>()
     )
    );
};

template<class ...Args>
using fom_jacobian_action_t = typename fom_jacobian_action<Args...>::type;

// -----------------------------------------------------------------------------

template<class MaskerType, class OperandType, class = void>
struct mask_action{
  using type = void;
};

template<class MaskerType, class OperandType>
struct mask_action<
  MaskerType, OperandType,
  mpl::enable_if_t<
    !std::is_void<
    decltype
    (std::declval<MaskerType const>().createApplyMaskResult
    (
    std::declval<OperandType const &>()
    )
    )
    >::value
  >
  >
{
  using type =
    decltype
    (std::declval<MaskerType const>().createApplyMaskResult
     (
      std::declval<OperandType const &>()
     )
    );
};

template<class ...Args>
using mask_action_t = typename mask_action<Args...>::type;

}}}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
