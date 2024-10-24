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

#ifndef PRESSIO_ROM_ROM_CONCEPTS_CXX20_HPP_
#define PRESSIO_ROM_ROM_CONCEPTS_CXX20_HPP_

#include <concepts>
#include "concepts_helpers.hpp"
#include "./impl/ode_has_const_discrete_residual_jacobian_action_method.hpp"

namespace pressio{ namespace rom{

template<class T>
concept ReducedState = ::pressio::is_vector_eigen<T>::value;

// ----------------------------------------------------------------------------
// SUBSPACES
// ----------------------------------------------------------------------------
template <class T>
concept VectorSubspace =
  requires(){ typename T::basis_matrix_type; }
  && std::copy_constructible<T>
  && std::copy_constructible<typename T::basis_matrix_type>
  && Traits<typename T::basis_matrix_type>::rank == 2
  && !std::assignable_from<T&, T>
  && !std::assignable_from<T&, T&>
  //
  && requires(const T & A)
  {
    { A.basis()         } -> std::same_as<const typename T::basis_matrix_type &>;
    { A.dimension()     } -> std::integral;
    { A.isColumnSpace() } -> std::convertible_to<bool>;
    { A.isRowSpace()    } -> std::convertible_to<bool>;
  };

template <class T>
concept PossiblyAffineTrialColumnSubspace =
  ReducedState<typename T::reduced_state_type>
  && std::copy_constructible< typename T::full_state_type>
  && std::copy_constructible< typename T::basis_matrix_type>
  && Traits<typename T::full_state_type>::rank   == 1
  && Traits<typename T::basis_matrix_type>::rank == 2
  /**/
  && std::copy_constructible<T>
  && !std::assignable_from<T&, T>
  && !std::assignable_from<T&, T&>
  //
  && requires(const T & A,
	      const typename T::reduced_state_type & reducedState,
	      typename T::full_state_type & fullState)
  {
    {A.createReducedState() } -> std::same_as<typename T::reduced_state_type>;
    {A.createFullState()    } -> std::same_as<typename T::full_state_type>;
    {A.createFullStateFromReducedState(reducedState) }
         -> std::same_as<typename T::full_state_type>;
  }
  //
  && requires(const T & A,
	      const typename T::reduced_state_type & reducedState,
	      typename T::full_state_type & fullState)
  {
    {A.mapFromReducedState(reducedState, fullState) } -> std::same_as<void>;
    {A.basisOfTranslatedSpace()} -> std::same_as<const typename T::basis_matrix_type &>;
    {A.translationVector()     } -> std::same_as<const typename T::full_state_type &>;
    {A.basis()        } -> std::same_as<const typename T::basis_matrix_type &>;
    {A.dimension()    } -> std::integral;
    {A.isColumnSpace()} -> std::convertible_to<bool>;
    {A.isRowSpace()   } -> std::convertible_to<bool>;
  };

template <class T>
concept PossiblyAffineRealValuedTrialColumnSubspace =
  PossiblyAffineTrialColumnSubspace<T>
  && std::floating_point< scalar_trait_t<typename T::reduced_state_type> >
  && std::floating_point< scalar_trait_t<typename T::full_state_type> >
  && std::floating_point< scalar_trait_t<typename T::basis_matrix_type> >;


// ----------------------------------------------------------------------------
// VARIOUS
// ----------------------------------------------------------------------------
template <class T, class MaskerType>
concept MaskableWith =
  requires(const T & operand, const MaskerType & masker)
  {
    { masker.createResultOfMaskActionOn(operand) } -> std::copy_constructible;
    { masker(operand,
	     std::declval<decltype(masker.createResultOfMaskActionOn(operand)) &>()
	     ) } -> std::same_as<void>;
  };

// ----------------------------------------------------------------------------
// SYSTEMS
// ----------------------------------------------------------------------------
template <class T, class JacobianActionOperandType>
concept SteadyFomWithJacobianAction =
     std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::residual_type>
  && requires(const T & A,
	      const typename T::state_type   & state,
	      const JacobianActionOperandType & operand)
  {
    { A.createResidual() } -> std::same_as<typename T::residual_type>;
    { A.createResultOfJacobianActionOn(operand) } -> std::copy_constructible;
  }
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::residual_type & residual,
	      const JacobianActionOperandType & operand,
	      std::optional<
	         impl::fom_jac_action_t<T, JacobianActionOperandType> *
	      > jActionOp)
  {
    { A.residualAndJacobianAction(state, residual,
				  operand, jActionOp) } -> std::same_as<void>;
  };

template <class T>
concept SemiDiscreteFom =
     std::regular<typename T::time_type>
  && std::totally_ordered<typename T::time_type>
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && requires(const T & A,
	      const typename T::state_type     & state,
	      const typename T::time_type      & evalTime,
	      typename T::rhs_type & rhs)
  {
    { A.createRhs() } -> std::same_as<typename T::rhs_type>;
    { A.rhs(state, evalTime, rhs) } -> std::same_as<void>;
  };

template <class T, class MassMatrixActionOperandType>
concept SemiDiscreteFomWithMassMatrixAction =
  SemiDiscreteFom<T>
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const MassMatrixActionOperandType & operand)
  {
    { A.createResultOfMassMatrixActionOn(operand) } -> std::copy_constructible;
  }
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const typename T::time_type     & evalTime,
	      const MassMatrixActionOperandType & operand,
	      impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType> & result)
  {
    { A.applyMassMatrix(state, operand, evalTime, result) } -> std::same_as<void>;
  };

template<class T, class JacobianActionOperandType>
concept SemiDiscreteFomWithJacobianAction =
  SemiDiscreteFom<T>
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const JacobianActionOperandType & operand)
  {
    { A.createResultOfJacobianActionOn(operand) } -> std::copy_constructible;
  }
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const typename T::time_type     & evalTime,
	      const JacobianActionOperandType & operand,
	      impl::fom_jac_action_t<T, JacobianActionOperandType> & result)
  {
    { A.applyJacobian(state, operand, evalTime, result) } -> std::same_as<void>;
  };

template<class T, class OperandType>
concept SemiDiscreteFomWithJacobianAndMassMatrixAction =
     SemiDiscreteFomWithJacobianAction<T, OperandType>
  && SemiDiscreteFomWithMassMatrixAction<T, OperandType>;


template<class T, int TotalNumStates, class JacobianActionOperandType>
concept FullyDiscreteSystemWithJacobianAction =
  (TotalNumStates == 2 || TotalNumStates == 3)
  && std::regular<typename T::time_type>
  && std::totally_ordered<typename T::time_type>
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::discrete_residual_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const JacobianActionOperandType & operand)
  {
    { A.createDiscreteTimeResidual() } -> std::same_as<typename T::discrete_residual_type>;
    { A.createResultOfDiscreteTimeJacobianActionOn(operand) } -> std::copy_constructible;
  }
  /*todo: fix syntax */
  && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
      T, TotalNumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::time_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      JacobianActionOperandType,
      impl::fully_discrete_fom_jac_action_t<T, JacobianActionOperandType>
    >::value;

// ----------------------------------------------------------------------------
// SYSTEMS REFINEMENTS FOR REAL VALUED
// ----------------------------------------------------------------------------
template <class T, class JacobianActionOperandType>
concept RealValuedSteadyFomWithJacobianAction =
  SteadyFomWithJacobianAction<T, JacobianActionOperandType>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::residual_type> >
  && std::floating_point<
       scalar_trait_t< impl::fom_jac_action_t<T, JacobianActionOperandType> >
     >;

template <class T>
concept RealValuedSemiDiscreteFom =
  SemiDiscreteFom<T>
  && std::floating_point< typename T::time_type>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >;

template <class T, class MassMatrixActionOperandType>
concept RealValuedSemiDiscreteFomWithMassMatrixAction =
  RealValuedSemiDiscreteFom<T>
  && SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>
  && std::floating_point<
    scalar_trait_t< impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType> >
  >;

template <class T, class JacobianActionOperandType>
concept RealValuedSemiDiscreteFomWithJacobianAction =
  RealValuedSemiDiscreteFom<T>
  && SemiDiscreteFomWithJacobianAction<T, JacobianActionOperandType>
  && std::floating_point<
       scalar_trait_t< impl::fom_jac_action_t<T, JacobianActionOperandType> >
     >;

template<class T, class OperandType>
concept RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction =
  RealValuedSemiDiscreteFomWithJacobianAction<T, OperandType>
  && RealValuedSemiDiscreteFomWithMassMatrixAction<T, OperandType>;

template<class T, int TotalNumStates, class JacobianActionOperandType>
concept RealValuedFullyDiscreteSystemWithJacobianAction =
  FullyDiscreteSystemWithJacobianAction<T, TotalNumStates, JacobianActionOperandType>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::discrete_residual_type> >
  && std::floating_point<
       scalar_trait_t< impl::fully_discrete_fom_jac_action_t<T, JacobianActionOperandType> >
     >;

}} // end namespace pressio::rom

#endif  // PRESSIO_ROM_ROM_CONCEPTS_CXX20_HPP_
