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

#include <typeinfo>

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

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
  using sc_t		 = typename T::sc_t;
  using solver_tag	 = solvertag;
  using this_t		 = Solver<solver_tag, T>;
  using state_t		 = typename T::state_t;
  using iterative_base_t = IterativeBase<this_t>;
  friend iterative_base_t;
  using typename iterative_base_t::iteration_t;

private:
  const sc_t defaultTol_  = static_cast<sc_t>(0.000001);

  iteration_t iStep_ = {};
  iteration_t jacobianUpdateFreq_ = 1;

  //0: CorrectionAbsoluteNorm
  //1: CorrectionRelativeNorm
  //2: ResidualAbsoluteNorm
  //3: ResidualRelativeNorm
  //4: GradientAbsoluteNorm
  //5: GradientRelativeNorm
  std::array<sc_t, 6> norms_ = {};

  //0: tol for CorrectionAbsoluteNorm
  //1: tol for CorrectionRelativeNorm
  //2: tol for ResidualAbsoluteNorm
  //3: tol for ResidualRelativeNorm
  //4: tol for GradientAbsoluteNorm
  //5: tol for GradientRelativeNorm
  std::array<sc_t, 6> tolerances_ = {};

  // updating criterion enum
  update updatingE_ = update::standard;
  std::shared_ptr<impl::BaseUpdater> updater_ = nullptr;

  // stopping creterion enum
  stop stoppingE_   = stop::whenCorrectionAbsoluteNormBelowTolerance;

public:
  Solver() = delete;
  Solver(Solver const &) = default;
  Solver & operator=(Solver const &) = delete;
  Solver(Solver &&) = default;
  Solver & operator=(Solver &&) = delete;
  ~Solver() = default;

  template <typename system_t, typename state_t, typename ...Args>
  Solver(const system_t & system,
	 const state_t & state,
	 stop stoppingE,
	 update updatingE,
	 Args &&... args)
    : T(system, state, std::forward<Args>(args)...)
    , updatingE_(updatingE)
    , stoppingE_(stoppingE)
  {
    tolerances_.fill(defaultTol_);
  }

  template <typename system_t, typename state_t, typename ...Args>
  Solver(const system_t & system,
	 const state_t & state,
	 Args &&... args)
    : Solver(system, state,
	     stop::whenCorrectionAbsoluteNormBelowTolerance,
	     update::standard,
	     std::forward<Args>(args)...)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /*  here we use this trick just to simplify code for
      pressio4py so that users can pass a ROM problem directly.
      But this is only supposed to be enabled when
      doing bindings for pressio4py */
  template <
    typename rom_problem_t, typename state_t, typename ...Args,
    mpl::enable_if_t<has_system_ref_method<rom_problem_t>::value, int> = 0
    >
  Solver(rom_problem_t & problem,
	 const state_t & state,
	 Args && ...args)
    : Solver(problem.systemRef(),state,std::forward<Args>(args)...)
  {}

  template <
    typename rom_problem_t, typename state_t, typename ...Args,
    mpl::enable_if_t<has_stepper_ref_method<rom_problem_t>::value, int> = 0
    >
  Solver(rom_problem_t & problem,
	 const state_t & state,
	 Args && ...args)
    : Solver(problem.stepperRef(),state,std::forward<Args>(args)...)
  {}
#endif

public:
  void setSystemJacobianUpdateFreq(std::size_t newFreq){
    jacobianUpdateFreq_ = newFreq;
  }

  iteration_t numIterationsExecuted() const {
    return iStep_;
  }

  // *****************************************
  // *** set or query updatating criterion ***
  void setUpdatingCriterion(update value){
    updatingE_ = value;
    // set null to indicate it needs to be constructed
    updater_ = nullptr;
  }

  update updatingCriterion() const{
    return updatingE_;
  }

  // *****************************************
  // *** set or query stopping criterion ***
  void setStoppingCriterion(stop value){
    stoppingE_ = value;
  }
  stop stoppingCriterion() const{
    return stoppingE_;
  }

  // *****************************************
  // *** set or query tolerances ***

  // this is used to set a single tol for all
  void setTolerance(sc_t tolerance){ tolerances_.fill(std::move(tolerance)); }

  // finer-grained methods for tolerances
  void setCorrectionAbsoluteTolerance(sc_t value){ tolerances_[0] = std::move(value); }
  void setCorrectionRelativeTolerance(sc_t value){ tolerances_[1] = std::move(value); }
  void setResidualAbsoluteTolerance(sc_t value)	 { tolerances_[2] = std::move(value); }
  void setResidualRelativeTolerance(sc_t value)  { tolerances_[3] = std::move(value); }
  void setGradientAbsoluteTolerance(sc_t value)  { tolerances_[4] = std::move(value); }
  void setGradientRelativeTolerance(sc_t value)  { tolerances_[5] = std::move(value); }

  sc_t correctionAbsoluteTolerance()const { return tolerances_[0]; }
  sc_t correctionRelativeTolerance()const { return tolerances_[1]; }
  sc_t residualAbsoluteTolerance()const   { return tolerances_[2]; }
  sc_t residualRelativeTolerance()const   { return tolerances_[3]; }
  sc_t gradientAbsoluteTolerance()const   { return tolerances_[4]; }
  sc_t gradientRelativeTolerance()const   { return tolerances_[5]; }

  template<typename system_t>
  void solve(const system_t & system, state_t & state)
  {
    this->solveImpl(system, state);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // this overload is needed when calling the solver from Python directly,
  // For example, this happens when doing steady LSPG since
  // the state vector is owned by the Python code
  template<typename system_t, typename _state_t = state_t>
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_array_pybind<_state_t>::value
    >
  solve(const system_t & system, _state_t & state)
  {
    // here we want to view the state since we want to modify its data,
    // which is numpy array owned by the user inside their Python code.
    // upon exit of this function, the original state is changed.
    ::pressio::containers::Tensor<1, _state_t> stateView(state, ::pressio::view());
    this->solveImpl(system, stateView);
  }
#endif

private:
  template<typename system_t, typename state_t>
  void solveImpl(const system_t & system, state_t & state)
  {
    PRESSIOLOG_INFO("nonlinsolver: solve");

    if (!updater_){
      PRESSIOLOG_DEBUG("nonlinsolver: creating updater");
      updater_ = createUpdater<solvertag>(state, updatingE_);
    }
    // after construction, it should NOT be null
    assert(updater_);

    sc_t residualNorm0 = {};
    sc_t correctionNorm0 = {};
    sc_t gradientNorm0 = {};
    bool recomputeSystemJacobian = true;

    iStep_ = 0;
    while (++iStep_ <= iterative_base_t::maxIters_)
    {
      recomputeSystemJacobian =
	(iStep_ == 1) ? true : ((iStep_ % jacobianUpdateFreq_) == 0);

      // 1.
      try{
	T::computeCorrection(system, state, recomputeSystemJacobian);
      }
      catch (::pressio::eh::residual_evaluation_failure_unrecoverable const &e)
      {
	PRESSIOLOG_CRITICAL("nonlinsolver: failure");
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
	  (iStep_,
	   norms_[0], norms_[1], norms_[2],
	   norms_[3], norms_[4], norms_[5]);
      }
      else{
	impl::printMetrics
	  (iStep_,
	   norms_[0], norms_[1], norms_[2], norms_[3]);
      }

      // 4.
      if (stopLoop(iStep_)){
	PRESSIOLOG_DEBUG("nonlinsolver: stopping");
	break;
      }

      // 5.
      applyUpdater(system, state, *this, updatingE_, updater_);
    }

    // when we are done with a solver, reset params to default
    T::resetForNewCall();
    updater_->resetForNewCall();
  }

  // stopping check
  bool stopLoop(const iteration_t & iStep) const
  {
    switch (stoppingE_)
      {
      case stop::afterMaxIters:
    	return iStep == iterative_base_t::maxIters_;

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

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
