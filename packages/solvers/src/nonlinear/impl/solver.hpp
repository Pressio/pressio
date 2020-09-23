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

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename solver_tag, typename T, typename sc_t>
class Solver
  : public T,
    public IterativeBase< Solver<solver_tag, T, sc_t>, sc_t >
{
  using this_t		 = Solver<solver_tag, T, sc_t>;

  using state_t		 = typename T::state_t;
  using iterative_base_t = IterativeBase<this_t, sc_t>;
  friend iterative_base_t;
  using typename iterative_base_t::iteration_t;

private:
  iteration_t jacobianUpdateFreq_ = 1;
  iteration_t iStep_ = {};
  std::array<sc_t, 6> norms_;

  // default stopping creterion
  stop stopping_   = stop::whenCorrectionAbsoluteNormBelowTolerance;

  // updating criterion
  update updating_ = update::standard;
  using upd_base_t = impl::BaseUpdater;
  using upd_def_t  = impl::DefaultUpdater;
  std::unique_ptr<upd_base_t> updater_ = nullptr;

public:
  Solver() = delete;

  template <typename system_t, typename state_t, typename ...Args>
  Solver(const system_t & system,
  	 const state_t & state,
	 stop stoppingE,
	 update updatingE,
  	 Args &&... args)
    : T(system, state, std::forward<Args>(args)...),
      stopping_(stoppingE),
      updating_(updatingE)
  {}

  template <typename system_t, typename state_t, typename ...Args>
  Solver(const system_t & system,
	 const state_t & state,
  	 Args &&... args)
    : Solver(system, state,
	     stop::whenCorrectionAbsoluteNormBelowTolerance,
	     update::standard,
	     std::forward<Args>(args)...)
  {}

  // copy constr and assign
  Solver(Solver const &) = default;
  Solver & operator=(Solver const &) = default;

  // move constr and assign
  Solver(Solver &&) = default;
  Solver & operator=(Solver &&) = default;

  // destr
  ~Solver() = default;

public:
  void setStoppingCriterion(stop value){
    stopping_ = value;
  }

  stop getStoppingCriterion() const{
    return stopping_;
  }

  void setUpdatingCriterion(update value){
    updating_ = value;
    // set null to indicate it needs to be constructed
    updater_ = nullptr;
  }

  update getUpdatingCriterion() const{
    return updating_;
  }

  void setSystemJacobianUpdateFreq(std::size_t newFreq){
    jacobianUpdateFreq_ = newFreq;
  }

public:

  template<typename system_t, typename _state_t = state_t>
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  mpl::enable_if_t<
  !::pressio::containers::predicates::is_array_pybind<_state_t>::value
  >
#else
  void
#endif
  solve(const system_t & sys, _state_t & state)
  {
    this->solveImpl(sys, state);
  }


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_array_pybind<state_t>::value
    >
  solve(const system_t & sys, state_t & state)
  {
    // here we want to view the state since we want to modify its data,
    // which is numpy array owned by the user inside their Python code.
    // upon exit of this function, the original state is changed.
    ::pressio::containers::Vector<state_t> stateView(state, ::pressio::view());
    this->solveImpl(sys, stateView);
  }
#endif

private:
  template<typename system_t, typename state_t>
  void solveImpl(const system_t & sys, state_t & state)
  {
    if (!updater_) constructUpdater(state);
    // after we construct it, it should not be nullt
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
	T::computeCorrection(sys, state, recomputeSystemJacobian);
      }
      catch (::pressio::eh::residual_evaluation_failure_unrecoverable const &e){
	throw ::pressio::eh::nonlinear_solve_failure();
      }

      // 2.
      const auto correctionNorm = T::correctionNormCurrentCorrectionStep();
      const auto residualNorm	= T::residualNormCurrentCorrectionStep();
      if (iStep_==1) {
	residualNorm0   = residualNorm;
	correctionNorm0 = correctionNorm;
      }

      norms_[0] = correctionNorm;
      norms_[1] = correctionNorm/correctionNorm0;
      norms_[2] = residualNorm;
      norms_[3] = residualNorm/residualNorm0;

      if (T::hasGradientComputation()){
	const auto gradientNorm	= T::gradientNormCurrentCorrectionStep();
	if (iStep_==1) gradientNorm0 = gradientNorm;

	norms_[4] = gradientNorm;
	norms_[5] = gradientNorm/gradientNorm0;
      }

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      if (T::hasGradientComputation()){
	impl::printNonlinearLeastSquaresDefaultMetrics
	  (iStep_, norms_[0], norms_[1], norms_[2],
	   norms_[3], norms_[4], norms_[5]);
      }
      else{
	impl::printNonlinearLeastSquaresDefaultMetrics
	  (iStep_, norms_[0], norms_[1], norms_[2], norms_[3]);
      }
#endif

      // 4.
      if (stopLoop(iStep_)) break;

      // 5.
      executeUpdater(sys, state);
    }

    // when we are done with a solver, reset params to default
    T::resetForNewCall();
    updater_->resetForNewCall();
  }

  iteration_t getNumIterationsExecutedImpl() const {
    return iStep_;
  }

  /*
    for GaussNewton, we can only have: standard, armijo update
  */
  template<typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<
    std::is_same<_solver_tag, GaussNewton>::value or
    std::is_same<_solver_tag, GaussNewtonQR>::value
    >
  constructUpdater(const state_t & state)
  {
    switch (updating_)
    {
    case update::standard:{
      using ut = impl::DefaultUpdater;
      updater_.reset(new ut());
      break;
    }

    case update::armijo:{
      using ut = impl::ArmijoUpdater<state_t>;
      updater_.reset(new ut(state));
      break;
    }

    default:
      throw std::runtime_error("Invalid update enum for GaussNewton");
    }
  }

  template<typename system_t, typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<
    std::is_same<_solver_tag, GaussNewton>::value or
    std::is_same<_solver_tag, GaussNewtonQR>::value
    >
  executeUpdater(const system_t & sys, state_t & state)
  {
    switch (updating_)
    {
      case update::standard:{
  	using ut = impl::DefaultUpdater;
  	updater_->template updateState<ut>(sys, state, *this);
  	break;
      }

      case update::armijo:{
  	using ut = impl::ArmijoUpdater<state_t>;
  	updater_->template updateState<ut>(sys, state, *this);
  	break;
      }

      default:
	throw std::runtime_error("Invalid update enum for GaussNewton");
    }
  }


  /*
    Levenberg-Mqrquardt, we can only have: standard, LM1, LM2
  */
  template<typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<std::is_same<_solver_tag, LevenbergMarquardt>::value>
  constructUpdater(const state_t & state)
  {
    switch (updating_)
    {
      case update::standard:{
	using ut = impl::DefaultUpdater;
	updater_.reset(new ut());
	break;
      }

      case update::LMSchedule1:{
	using ut = impl::LMSchedule1Updater<state_t>;
	updater_.reset(new ut(state));
	break;
      }

      case update::LMSchedule2:{
	using ut = impl::LMSchedule2Updater<state_t>;
	updater_.reset(new ut(state));
	break;
      }

      default:
	throw std::runtime_error("Invalid update enum for LevenbergMarquardt");
    }
  }

  template<typename system_t, typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<std::is_same<_solver_tag, LevenbergMarquardt>::value>
  executeUpdater(const system_t & sys, state_t & state)
  {
    switch (updating_)
      {
      case update::standard:{
	using ut = impl::DefaultUpdater;
	updater_->template updateState<ut>(sys, state, *this);
	break;
      }

      case update::LMSchedule1:{
	using ut = impl::LMSchedule1Updater<state_t>;
	updater_->template updateState<ut>(sys, state, *this);
	break;
      }

      case update::LMSchedule2:{
	using ut = impl::LMSchedule2Updater<state_t>;
	updater_->template updateState<ut>(sys, state, *this);
	break;
      }

      default:
	throw std::runtime_error("Invalid update enum for LevenbergMarquardt");
      }
  }


  /*
    Newton-Raphson, we can only have: standard
  */
  template<typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<std::is_same<_solver_tag, NewtonRaphson>::value>
  constructUpdater(const state_t & state)
  {
    switch (updating_)
      {
      case update::standard:{
  	using ut = impl::DefaultUpdater;
  	updater_.reset(new ut());
  	break;
      }

      default:
  	throw std::runtime_error("Invalid update enum for NewtonRaphson");
      }
  }

  template<typename system_t, typename state_t, typename _solver_tag = solver_tag>
  mpl::enable_if_t<std::is_same<_solver_tag, NewtonRaphson>::value>
  executeUpdater(const system_t & sys, state_t & state)
  {
    switch (updating_)
      {
      case update::standard:{
  	using ut = impl::DefaultUpdater;
  	updater_->template updateState<ut>(sys, state, *this);
  	break;
      }
      default:
  	throw std::runtime_error("Invalid update enum for NewtonRaphson");
      }
  }


  // stopping check
  bool stopLoop(const iteration_t & iStep) const
  {
    switch (stopping_)
      {
      case stop::afterMaxIters:
    	return iStep == iterative_base_t::maxIters_;

      case stop::whenCorrectionAbsoluteNormBelowTolerance:
    	return norms_[0] < iterative_base_t::tolerance_;
      case stop::whenCorrectionRelativeNormBelowTolerance:
    	return norms_[1] < iterative_base_t::tolerance_;

      case stop::whenResidualAbsoluteNormBelowTolerance:
    	return norms_[2] < iterative_base_t::tolerance_;
      case stop::whenResidualRelativeNormBelowTolerance:
    	return norms_[3] < iterative_base_t::tolerance_;

      case stop::whenGradientAbsoluteNormBelowTolerance:
    	return norms_[4] < iterative_base_t::tolerance_;
      case stop::whenGradientRelativeNormBelowTolerance:
    	return norms_[5] < iterative_base_t::tolerance_;

      default:
    	return false;
      };
  }

};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
