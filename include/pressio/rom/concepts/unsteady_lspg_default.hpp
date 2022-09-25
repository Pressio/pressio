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

namespace pressio{ namespace rom{ namespace lspg{

#ifdef PRESSIO_ENABLE_CXX20

template <class T, class SubspaceType>
concept DefaultDiscreteTimeAssemblyWith =
      SemiDiscreteFomWithJacobianAction<T, SubspaceType>
   && PossiblyAffineTrialColumnSubspace<SubspaceType>
   && requires(const SubspaceType & subspace,
	       const T & A,
	       const typename Traits<typename T::state_type>::scalar_type & scalarCoeff,
               const typename SubspaceType::basis_matrix_type & basisMatrix)
  {

    ::pressio::ops::deep_copy(std::declval<typename T::state_type &>(),
			      std::declval<typename T::state_type const &>());
    ::pressio::ops::set_zero(std::declval<typename T::state_type &>());

    // needed by ode rhs manager
    ::pressio::ops::set_zero(std::declval<typename T::state_type &>());

    // needed by lspg impl
    ::pressio::ops::update(std::declval<decltype(A.createApplyJacobianResult(basisMatrix)) &>(),
			   scalarCoeff, basisMatrix, scalarCoeff);

    ::pressio::ops::update(std::declval<typename T::right_hand_side_type &>(), scalarCoeff,
			   std::declval<typename T::state_type const &>(), scalarCoeff,
			   std::declval<typename T::state_type const &>(), scalarCoeff);

    ::pressio::ops::update(std::declval<typename T::right_hand_side_type &>(), scalarCoeff,
			   std::declval<typename T::state_type const &>(), scalarCoeff,
			   std::declval<typename T::state_type const &>(), scalarCoeff,
			   std::declval<typename T::state_type const &>(), scalarCoeff);
  };



#else

template<class T, class SubspaceType, class enable = void>
struct DefaultDiscreteTimeAssemblyWith : std::false_type{};

template<class T, class SubspaceType>
struct DefaultDiscreteTimeAssemblyWith<
  T, SubspaceType,
  mpl::enable_if_t<
       SemiDiscreteFomWithJacobianAction<T, SubspaceType>::value
    && PossiblyAffineTrialColumnSubspace<SubspaceType>::value
   >
  > : std::true_type{};

#endif

}}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
