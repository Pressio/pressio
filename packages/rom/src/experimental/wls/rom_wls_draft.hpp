/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_system.hpp
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

#ifndef ROM_WLS_PROBLEM_RESIDUAL_API_HPP_
#define ROM_WLS_PROBLEM_RESIDUAL_API_HPP_

namespace pressio{ namespace rom{ namespace experimental{

// container for the fom states for WLS
template <std::size_t n, typename fom_state_type, typename reconstuctor_type>
class WlsFomStatesContainer{

  // put here usual things and overload operator [] so we can access the fom state at a given index
  // where the index typically is the step number
public:

  static constexpr std::size_t n_ = n;

  WlsFomStatesContainer() = delete;
  ~WlsFomStatesContainer() = default;

  template <
    typename _fom_state_type = fom_state_type,
    ::pressio::mpl::enable_if_t<
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  WlsFomStatesContainer(const _fom_state_type & fomStateIn,
			const reconstuctor_type & fomStateReconstr)
    : fomStates_{fomStateIn},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

public:
  fom_state_type & operator[](std::size_t index) {
    return fomStates_[index];
  }

  fom_state_type const & operator[](std::size_t index) const{
    return fomStates_[index];
  }

  template <typename rom_state_t>
  void reconstructFomStateAt(const rom_state_t & romY, std::size_t index) const
  {
    fomStateReconstrObj_(romY[index], fomStates_[index]);
  }

private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    for (auto i=0; i<n; i++)
      ::pressio::containers::ops::set_zero(fomStates_[i]);
  }

private:
  mutable std::array<fom_state_type, n> fomStates_;
  const reconstuctor_type & fomStateReconstrObj_  = {};
};




/* the WlsSystem is the class that allows to interface with the solvers */
template<
  typename fom_type,
  typename residual_policy_type,
  typename jacobian_policy_type,
  typename wls_state_type,
  typename wls_residual_type,
  typename wls_jacobian_type
  >
class WlsSystem{

  const fom_type & fomObj_;
  const residual_policy_type & residualPolicyObj_;
  const jacobian_policy_type & jacobianPolicyObj_;

public:
  // these need to be public because are detected by solver
  using scalar_type	= typename fom_type::scalar_type;
  using state_type	= wls_state_type;
  using residual_type	= wls_residual_type;
  using jacobian_type	= wls_jacobian_type;

public:
  WlsSystem() = delete;
  ~WlsSystem() = default;

  WlsSystem(const fom_type & fomObject,
	    const residual_policy_type & resPolicyObj,
	    const jacobian_policy_type & jacPolicyObj)
    : fomObj_(fomObject),
      residualPolicyObj_(resPolicyObj),
      jacobianPolicyObj_(jacPolicyObj){}

public:
  void residual(const wls_state_type & wlsState, wls_residual_type & R) const{
    (this->residualPolicyObj_).template operator()(wlsState, fomObj_, R);
  }
]]
  // void jacobian(const wls_state_type & wlsState, wls_jacobian_type & J) const{
  //   // call somehow jacobian policy
  // }

  wls_residual_type residual(const wls_state_type & wlsState) const{
    return (this->residualPolicyObj_).template operator()(wlsState, fomObj_);
  }

  // wls_jacobian_type jacobian(const wls_state_type & wlsState) const{
  //   return // call somehow jacobian policy
  // }
};//end class





// n is the window size
template <std::size_t n, typename fom_type, typename wls_state_type, typename ... Args>
struct DefaultWlsTypeGeneratorResidualApi{

  // native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename fom_t::scalar_type;
  using fom_native_state_t	= typename fom_t::state_type;
  using fom_native_residual_t	= typename fom_t::residual_type;

  // fom_state_t: this is a wrapper of the fom native state, and using a vector wrapper
  // seems like a good idea since typically the fom uses a vector for the state
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;

  // fom_residual_t: this is a wrapper of the fom native residual, and using a vector wrapper
  // seems like a good idea since typically the fom uses a vector for the residual
  using fom_residual_t		= ::pressio::containers::Vector<fom_native_residual_t>;

  // wls_state_t is the type holding the wls state, so maybe this
  // is a multi-vector or some custom type that we make. On object of this
  // is supposed to own the generalized coordinates at the n steps of the target window.
  // the wls_state_type is passed by the user so it is defined in the main file.
  using wls_state_t		= wls_state_type;

  // wls_residual_t is the type holding the wls residual, so maybe this should be made
  // a multi-vector of fom_residual_t or some custom type that we make.
  // For this DefaultWlsTypeGenerator, we can make this a multi-vector.
  // If we need something else, we can create a different WlsTypeGenerator.
  using wls_residual_t		= ::pressio::containers::MultiVector<fom_residual_t>;

  // wls_matrix_t: an object of this type should hold somehow the wls matrix, so basically the large
  // J*phi that stems from the wls formulation. Since the WLS matrix has a block structure with dense blocks,
  // a basic version could be one where each block is a multivector (as it is done now for lspg) and
  // we use a std::list of multi-vectors to hold the full matrix.
  using wls_matrix_t		= /* */;

  // decoder types are easy, see the LspgUnsteady classes
  using decoder_t		= /* get it from args */;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<fom_state_t, decoder_t>;

  // fom_states_data: is used to store the fom states, use a custom class for this (see top of this file)
  using fom_states_container_t	= WlsFomStatesContainer<n, fom_state_t, fom_state_reconstr_t>;

  // here we are relying on the residual API, so we expect the user to tell us some details
  // on the stepper scheme used by the fom
  // --- find the total number of states needed ---
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::impl::IsStepperTotalNumStatesSetter, Args...>;
  using tot_n_setter = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert( !std::is_void<tot_n_setter>::value, "...");
  // numAuxStates is the number of auxiliary states needed by the fom object to
  // compute time discrete operators. This can be extracted from the args... like we do for lspg_unsteady
  static constexpr std::size_t auxStates = tot_n_setter::value - 1;

  // policy to compute the WLS residual: here we use the most basic one that computes things sequentially.
  // For other types of policies, we can create other type generators.
  using wls_residual_policy_t	= ::pressio::rom::experimental::WlsSequentialResidualPolicyForResidualApi<
    auxStates, wls_residual_t, fom_states_container_t>;

  // // policy to compute the WLS jacobian, do something similar for the residual one
  //using wls_jacobian_policy_t	= /* */;

};//end class



template<
  template <std::size_t n, typename, typename, typename ...> class wls_type,
  std::size_t n,
  typename fom_type,
  typename wls_state_type,
  typename ... Args>
class WlsProblemGeneratorResidualApi
{

public:
  using wls_problem_t = wls_type<fom_type, wls_state_type, Args...>;

  using fom_t			= typename wls_problem_t::fom_t;
  using scalar_t		= typename wls_problem_t::scalar_t;
  using fom_native_state_t	= typename wls_problem_t::fom_native_state_t;
  using fom_native_residual_t	= typename wls_problem_t::fom_native_residual_t;

  using fom_state_t		= typename wls_problem_t::fom_state_t;
  using fom_residual_t		= typename wls_problem_t::fom_residual_t;

  using wls_state_t		= typename wls_problem_t::wls_state_t;
  using wls_residual_t		= typename wls_problem_t::wls_residual_t;
  using wls_matrix_t		= typename wls_problem_t::wls_matrix_t;

  using decoder_t		= typename wls_problem_t::decoder_t;
  using fom_state_reconstr_t	= typename wls_problem_t::fom_state_reconstr_t;
  using fom_states_container_t  = typename wls_problem_t::fom_states_container_t;

  using wls_residual_policy_t	= typename wls_problem_t::wls_residual_policy_t;
  using wls_jacobian_policy_t	= typename wls_problem_t::wls_jacobian_policy_t;
  using wls_system_t		= typename wls_problem_t::wls_system_t;

  wls_residual_policy_t		residualPolicy_;
  wls_jacobian_policy_t		jacobianPolicy_;
  wls_system_t			systemObj_;

  WlsProblemGeneratorResidualApi(const fom_t	 & fomObj,
				 const fom_native_state_t & fomNativeStateReference,
				 decoder_t	 & decoder,
				 wls_state_t	 & wlsInitialState,
				 scalar_t	t0)
    : /* construct things that we need, like the policies and the system
         basically here we should put what Eric's already started doing. 
	residualPolicy(...),
	jacobianPolicy(..),
	systemObject(fomObj, residualPolicy, jacobianPolicy)
      */
      {}

};//end class





int main()
{
  //
  // all the usual stuff needed
  //

  constexpr std::size_t numStepsInWindow = 5;

  // how many states are needed by the stencil used by the user in the fom
  // for bdf1, this should be 2 because it uses y_n and y_n-1
  // for bdf2, this should be 3 because it uses y_n and y_n-1, y_n-2
  using stepper_stencil = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

  // the wls problem
  /* Solver Types
	1.) Regular Gauss--Newton
    a.) storage for N x Ns residuals (and associated routines)
    b.) Storage for a sparse (N Ns) x (K Ns ) Jacobian  
  */
  using wls_problem	 = pressio::rom::WlsProb lemGeneratorResidualApi<
    DefaultWlsTypeGeneratorResidualApi, numStepsInWindow,
    fom_t, wls_state_t, decoder_t, stepper_stencil, scalar_t>;

  wls_problem wlsProblem(appobj, /* whatever goes here */ );

}

}}} // end namespace pressio::rom::experimental
#endif
