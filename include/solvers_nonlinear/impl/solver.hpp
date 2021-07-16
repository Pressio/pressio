/*
//@HEADER
// ************************************************************************
//
// solver.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_

#include "solvers_iterative_base.hpp"
#include <typeinfo>

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename T, typename = void>
struct has_system_ref_method : std::false_type{};

template <typename T>
struct has_system_ref_method<
  T,
  mpl::enable_if_t<
    !std::is_void<
      decltype(std::declval<T>().systemRef())
      >::value
    >
  > : std::true_type{};


template <typename T, typename = void>
struct has_stepper_ref_method : std::false_type{};

template <typename T>
struct has_stepper_ref_method<
  T,
  mpl::enable_if_t<
    !std::is_void<
      decltype(std::declval<T>().stepperRef())
      >::value
    >
  > : std::true_type{};
#endif
// -------------------------------------------


template<typename solvertag, typename T>
class Solver
  : public T,
    public IterativeBase<Solver<solvertag, T>>
{
public:
  using sc_type		 = typename T::scalar_type;
  using solver_tag	 = solvertag;
  using this_type		 = Solver<solver_tag, T>;
  using state_type		 = typename T::state_type;
  using iterative_base_type = IterativeBase<this_type>;
  friend iterative_base_type;
  using typename iterative_base_type::iteration_t;

private:
  const sc_type defaultTol_  = static_cast<sc_type>(0.000001);

  iteration_t iStep_ = {};
  iteration_t jacobianUpdateFreq_ = 1;

  //0: CorrectionAbsoluteNorm
  //1: CorrectionRelativeNorm
  //2: ResidualAbsoluteNorm
  //3: ResidualRelativeNorm
  //4: GradientAbsoluteNorm
  //5: GradientRelativeNorm
  std::array<sc_type, 6> norms_ = {};

  //0: tol for CorrectionAbsoluteNorm
  //1: tol for CorrectionRelativeNorm
  //2: tol for ResidualAbsoluteNorm
  //3: tol for ResidualRelativeNorm
  //4: tol for GradientAbsoluteNorm
  //5: tol for GradientRelativeNorm
  std::array<sc_type, 6> tolerances_ = {};

  // updating criterion enum
  update updatingE_ = update::standard;
  std::unique_ptr<impl::BaseUpdater> updater_ = nullptr;

  // stopping creterion enum
  stop stoppingE_   = stop::whenCorrectionAbsoluteNormBelowTolerance;

  // size of system object is used to check if system changes
  // this might not be a great solution, but we will work on this
  std::size_t sizeOfSystemObj_ = {};

  // to monitor state during nonlinear solver
  std::unique_ptr<impl::BaseObserver> observer_ = nullptr;

  // by default, printing metrics is done with
  // column strings with the labels, but if one
  // wants to print only the numbers stripping everything else,
  // one can turn this on
  bool printStrippedMetrics_ = false;

public:
  Solver() = delete;
  Solver(Solver const &) = default;
  Solver & operator=(Solver const &) = delete;
  Solver(Solver &&) = default;
  Solver & operator=(Solver &&) = delete;
  ~Solver() = default;

  template <typename system_t, typename state_type, typename ...Args>
  Solver(const system_t & system,
	 const state_type & state,
	 stop stoppingE,
	 update updatingE,
	 Args &&... args)
    : T(system, state, std::forward<Args>(args)...)
    , updatingE_(updatingE)
    , stoppingE_(stoppingE)
  {
    tolerances_.fill(defaultTol_);
  }

  template <typename system_t, typename state_type, typename ...Args>
  Solver(const system_t & system,
	 const state_type & state,
	 Args &&... args)
    : Solver(system, state,
	     stop::whenCorrectionAbsoluteNormBelowTolerance,
	     update::standard,
	     std::forward<Args>(args)...)
  {}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   /*  here we use this trick just to simplify code for
//       pressio4py so that users can pass a ROM problem directly.
//       But this is only supposed to be enabled when
//       doing bindings for pressio4py */
//   template <
//     typename rom_problem_t, typename state_type, typename ...Args,
//     mpl::enable_if_t<has_system_ref_method<rom_problem_t>::value, int> = 0
//     >
//   Solver(rom_problem_t & problem,
// 	 const state_type & state,
// 	 Args && ...args)
//     : Solver(problem.systemRef(),state,std::forward<Args>(args)...)
//   {}

//   template <
//     typename rom_problem_t, typename state_type, typename ...Args,
//     mpl::enable_if_t<has_stepper_ref_method<rom_problem_t>::value, int> = 0
//     >
//   Solver(rom_problem_t & problem,
// 	 const state_type & state,
// 	 Args && ...args)
//     : Solver(problem.stepperRef(),state,std::forward<Args>(args)...)
//   {}
// #endif

public:
  void printStrippedMetrics(){
    printStrippedMetrics_ = true;
  }

  void setSystemJacobianUpdateFreq(std::size_t newFreq){
    jacobianUpdateFreq_ = newFreq;
  }

  iteration_t numIterationsExecuted() const {
    return iStep_;
  }

  // *** set or query updatating criterion ***
  void setUpdatingCriterion(update value){
    updatingE_ = value;
    // set 0 to indicate it needs to be constructed
    sizeOfSystemObj_ = 0;
  }

  update updatingCriterion() const{
    return updatingE_;
  }

  // *** set or query stopping criterion ***
  void setStoppingCriterion(stop value){
    stoppingE_ = value;
  }
  stop stoppingCriterion() const{
    return stoppingE_;
  }

  template<class functor_t>
  void setObserver(const functor_t & obsF)
  {
    PRESSIOLOG_INFO("nonlinsolver: creating observer");

    using o_t = Observer<state_type, const functor_t &>;
    observer_ = pressio::utils::make_unique<o_t>(obsF);
    observer_->applyFnc_ = applyObserver<o_t>;
  }

  // *** set or query tolerances ***

  // this is used to set a single tol for all
  void setTolerance(sc_type tolerance){ tolerances_.fill(std::move(tolerance)); }

  // finer-grained methods for tolerances
  void setCorrectionAbsoluteTolerance(sc_type value){ tolerances_[0] = std::move(value); }
  void setCorrectionRelativeTolerance(sc_type value){ tolerances_[1] = std::move(value); }
  void setResidualAbsoluteTolerance(sc_type value)	 { tolerances_[2] = std::move(value); }
  void setResidualRelativeTolerance(sc_type value)  { tolerances_[3] = std::move(value); }
  void setGradientAbsoluteTolerance(sc_type value)  { tolerances_[4] = std::move(value); }
  void setGradientRelativeTolerance(sc_type value)  { tolerances_[5] = std::move(value); }

  sc_type correctionAbsoluteTolerance()const { return tolerances_[0]; }
  sc_type correctionRelativeTolerance()const { return tolerances_[1]; }
  sc_type residualAbsoluteTolerance()const   { return tolerances_[2]; }
  sc_type residualRelativeTolerance()const   { return tolerances_[3]; }
  sc_type gradientAbsoluteTolerance()const   { return tolerances_[4]; }
  sc_type gradientRelativeTolerance()const   { return tolerances_[5]; }


  // *** solve ***
  template<typename system_t>
  void solve(const system_t & system, state_type & state)
  {
    // before we solve, we check if we need to recreate the updater
    // for example, this is the case if the system object changes
    if (recreateUpdater(system)){
      PRESSIOLOG_INFO("nonlinsolver: create updater");
      updater_ = createUpdater<this_type, system_t>(state, updatingE_);
    }
    this->solveImpl(system, state, *updater_);
  }

#if not defined PRESSIO_ENABLE_TPL_PYBIND11
  template<class system_t, class custom_updater_t>
  void solve(const system_t & system,
	     state_type & state,
	     custom_updater_t && updater)
  {
    // here we have to create the updater each time becuase
    // we don't know if the updater changes or not

    PRESSIOLOG_INFO("nonlinsolver: custom updater");
    updatingE_ = ::pressio::nonlinearsolvers::update::custom;

    using u_t =
      mpl::conditional_t<
	std::is_lvalue_reference<custom_updater_t&&>::value,
      Updater<system_t, state_type, this_type, custom_updater_t&>,
      Updater<system_t, state_type, this_type, custom_updater_t>
      >;
    updater_ = pressio::utils::make_unique<u_t>(std::forward<custom_updater_t>(updater));
    updater_->applyFnc_ = applyUpdater<u_t>;
    updater_->resetFnc_ = resetUpdater<u_t>;
    this->solveImpl(system, state, *updater_);
  }
#endif

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   // this overload is needed when calling the solver from Python directly,
//   // For example, this happens when doing steady LSPG since
//   // the state vector is owned by the Python code
//   template<typename system_t, typename _state_type = state_type>
//   mpl::enable_if_t<
//     ::pressio::is_array_pybind<_state_type>::value
//     >
//   solve(const system_t & system, _state_type & state)
//   {
//     // here we want to view the state since we want to modify its data,
//     // which is numpy array owned by the user inside their Python code.
//     // upon exit of this function, the original state is changed.
//     ::pressio::containers::Tensor<1, _state_type> stateView(state, ::pressio::view());

//     if (recreateUpdater(system)){
//       PRESSIOLOG_INFO("nonlinsolver: create updater");
//       updater_ = createUpdater<this_type, system_t>(stateView, updatingE_);
//     }
//     this->solveImpl(system, stateView, *updater_);
//   }
// #endif

private:
  template<typename system_t>
  bool recreateUpdater(const system_t & system)
  {
    if( (sizeOfSystemObj_ != sizeof(system)) or (sizeOfSystemObj_ == 0)){
      sizeOfSystemObj_ = sizeof(system);
      return true;
    }
    else if(!updater_){
      sizeOfSystemObj_ = sizeof(system);
      return true;
    }
    else{
      return false;
    }
  }

  template<class system_t, class state_type, class updater_t>
  void solveImpl(const system_t & system,
		 state_type & state,
		 updater_t & updater)
  {
    PRESSIOLOG_INFO("nonlinsolver: solve");

    sc_type residualNorm0 = {};
    sc_type correctionNorm0 = {};
    sc_type gradientNorm0 = {};
    bool recomputeSystemJacobian = true;

    iStep_ = 0;
    while (++iStep_ <= iterative_base_type::maxIters_)
    {
      recomputeSystemJacobian =
	(iStep_ == 1) ? true : ((iStep_ % jacobianUpdateFreq_) == 0);

      if (observer_){
	(*observer_)(iStep_, state);
      }

      // 1.
      try{
	T::computeCorrection(system, state, recomputeSystemJacobian);
      }
      catch (::pressio::eh::residual_evaluation_failure_unrecoverable const &e)
      {
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::nonlinear_solve_failure();
      }
      catch (::pressio::eh::residual_has_nans const &e)
      {
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::nonlinear_solve_failure();
      }

      // 2.
      const auto correctionNorm = T::correctionNormCurrentCorrectionStep();
      const auto residualNorm	= T::residualNormCurrentCorrectionStep();
      if (iStep_==1) {
	     residualNorm0   = residualNorm;
	     correctionNorm0 = correctionNorm;
      }
      norms_[0] = std::move(correctionNorm);
      norms_[1] = norms_[0]/correctionNorm0;
      norms_[2] = std::move(residualNorm);
      norms_[3] = norms_[2]/residualNorm0;

      if (T::hasGradientComputation()){
	const auto gradientNorm	= T::gradientNormCurrentCorrectionStep();
	if (iStep_==1) gradientNorm0 = gradientNorm;

	norms_[4] = std::move(gradientNorm);
	norms_[5] = gradientNorm/gradientNorm0;
      }

      if (T::hasGradientComputation()){
	impl::printMetrics
	  (iStep_, printStrippedMetrics_,
	   norms_[0], norms_[1], norms_[2],
	   norms_[3], norms_[4], norms_[5]);
      }
      else{
	impl::printMetrics
	  (iStep_, printStrippedMetrics_,
	   norms_[0], norms_[1], norms_[2], norms_[3]);
      }

      // 4.
      if (stopLoop(iStep_)){
	PRESSIOLOG_DEBUG("nonlinsolver: stopping");
	break;
      }

      // 5.
      updater(system, state, *this);
    }

    // when we are done with a solver, reset params to default
    T::resetForNewCall();
    updater_->reset();
  }

  // stopping check
  bool stopLoop(const iteration_t & iStep) const
  {
    switch (stoppingE_)
      {
      case stop::afterMaxIters:
    	return iStep == iterative_base_type::maxIters_;

      case stop::whenCorrectionAbsoluteNormBelowTolerance:
    	return norms_[0] < tolerances_[0];
      case stop::whenCorrectionRelativeNormBelowTolerance:
    	return norms_[1] < tolerances_[1];

      case stop::whenResidualAbsoluteNormBelowTolerance:
    	return norms_[2] < tolerances_[2];
      case stop::whenResidualRelativeNormBelowTolerance:
    	return norms_[3] < tolerances_[3];

      case stop::whenGradientAbsoluteNormBelowTolerance:
    	return norms_[4] < tolerances_[4];
      case stop::whenGradientRelativeNormBelowTolerance:
    	return norms_[5] < tolerances_[5];

      default:
    	return false;
      };
  }

};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
