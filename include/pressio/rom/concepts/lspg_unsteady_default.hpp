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

#ifndef ROM_CONSTRAINTS_ROM_UNSTEADY_LSPG_DEFAULT_HPP_
#define ROM_CONSTRAINTS_ROM_UNSTEADY_LSPG_DEFAULT_HPP_

#include "helpers.hpp"

#ifdef PRESSIO_ENABLE_CXX20

// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{

template <class TrialSubspaceType, class FomSystemType>
concept ComposableIntoDefaultProblem =
      SemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
   && PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
   && requires(const FomSystemType & A,
	       const typename TrialSubspaceType::basis_matrix_type & basisMatrix)
  {

    ::pressio::ops::deep_copy(std::declval<typename FomSystemType::state_type &>(),
			      std::declval<typename FomSystemType::state_type const &>());

 //    ::pressio::ops::set_zero(std::declval<typename T::state_type &>());
 //    ::pressio::ops::set_zero(std::declval<typename T::right_hand_side_type &>());

 //    ::pressio::ops::update
 //       (std::declval<typename T::right_hand_side_type &>(),
	// std::declval<scalar_trait_t<typename T::right_hand_side_type> const &>(),
	// std::declval<typename T::state_type const &>(),
	// std::declval<scalar_trait_t<typename T::state_type> const &>(),
	// std::declval<typename T::state_type const &>(),
	// std::declval<scalar_trait_t<typename T::state_type> const &>());

 //    ::pressio::ops::update
 //       (std::declval<typename T::right_hand_side_type &>(),
	// std::declval<scalar_trait_t<typename T::right_hand_side_type> const &>(),
	// std::declval<typename T::state_type const &>(),
	// std::declval<scalar_trait_t<typename T::state_type> const &>(),
	// std::declval<typename T::state_type const &>(),
	// std::declval<scalar_trait_t<typename T::state_type> const &>(),
	// std::declval<typename T::state_type const &>(),
	// std::declval<scalar_trait_t<typename T::state_type> const &>());

 //    ::pressio::ops::update
 //       (std::declval<decltype(A.createApplyJacobianResult(basisMatrix)) &>(),
	// std::declval<scalar_trait_t<typename TrialSubspaceType::basis_matrix_type> const &>(),
	// basisMatrix,
	// std::declval<scalar_trait_t<typename TrialSubspaceType::basis_matrix_type> const &>());
  };

}}}} //end namespace pressio::rom::lspg::unsteady

// #else

// namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{

// template<class T, class SubspaceType, class enable = void>
// struct DefaultDiscreteTimeAssemblyWith : std::false_type{};

// template<class T, class SubspaceType>
// struct DefaultDiscreteTimeAssemblyWith<
//   T, SubspaceType,
//   mpl::enable_if_t<
//        SemiDiscreteFomWithJacobianAction<T, SubspaceType>::value
//     && PossiblyAffineTrialColumnSubspace<SubspaceType>::value
//    >
//   > : std::true_type{};

// }}}}
#endif

#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
