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

#include <array>
#include "solvers_iterative_base.hpp"
#include <typeinfo>

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template<typename TagType, typename T>
class Solver
  : public T,
    public IterativeBase<Solver<TagType, T>>
{
public:
  using state_type = typename T::state_type;
  using residual_norm_type = typename T::residual_norm_type;
  using gradient_norm_type = typename T::gradient_norm_type;
  using correction_norm_type = typename T::correction_norm_type;

  using solver_tag = TagType;
  using this_type = Solver<solver_tag, T>;
  using iterative_base_type = IterativeBase<this_type>;
  friend iterative_base_type;
  using typename iterative_base_type::iteration_type;

private:

  iteration_type iStep_ = {};
  iteration_type jacobianUpdateFreq_ = 1;

  // CorrectionAbsoluteNorm, CorrectionRelativeNorm, and tolerances
  std::array<correction_norm_type, 2> correctionNorms_ = {};
  std::array<correction_norm_type, 2> correctionTolerances_ = {};

  // ResidualAbsoluteNorm, ResidualRelativeNorm, and tolerances
  std::array<residual_norm_type, 2> residualNorms_ = {};
  std::array<residual_norm_type, 2> residualTolerances_ = {};

  // GradientAbsoluteNorm, GradientRelativeNorm, and tolerances
  std::array<gradient_norm_type, 2> gradientNorms_ = {};
  std::array<gradient_norm_type, 2> gradientTolerances_ = {};

  // updating criterion enum
  Update updatingE_ = Update::Standard;
  std::unique_ptr<impl::BaseUpdater> updater_ = nullptr;

  // stopping creterion enum
  Stop stoppingE_   = Stop::WhenCorrectionAbsoluteNormBelowTolerance;

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

  template <typename SystemType, typename ...Args>
  Solver(const SystemType & system,
	 Stop stoppingE,
	 Update updatingE,
	 Args &&... args)
    : T(system, std::forward<Args>(args)...)
    , updatingE_(updatingE)
    , stoppingE_(stoppingE)
  {
    fillTolerancesWithDefaultValues();
  }

  template <typename SystemType, typename ...Args>
  Solver(const SystemType & system,
	 Args &&... args)
    : Solver(system,
	     Stop::WhenCorrectionAbsoluteNormBelowTolerance,
	     Update::Standard,
	     std::forward<Args>(args)...)
  {}

public:
  void printStrippedMetrics(){
    printStrippedMetrics_ = true;
  }

  void setSystemJacobianUpdateFreq(std::size_t newFreq){
    jacobianUpdateFreq_ = newFreq;
  }

  iteration_type numIterationsExecuted() const {
    return iStep_;
  }

  // *** set or query updatating criterion ***
  void setUpdatingCriterion(Update value){
    updatingE_ = value;
    // set 0 to indicate it needs to be constructed
    sizeOfSystemObj_ = 0;
  }

  Update updatingCriterion() const{
    return updatingE_;
  }

  // *** set or query stopping criterion ***
  void setStoppingCriterion(Stop value){
    stoppingE_ = value;
  }
  Stop stoppingCriterion() const{
    return stoppingE_;
  }

  template<class FunctorType>
  void setObserver(const FunctorType & obsF)
  {
    PRESSIOLOG_INFO("nonlinsolver: creating observer");

    using o_t = Observer<state_type, const FunctorType &>;
    observer_ = pressio::utils::make_unique<o_t>(obsF);
    observer_->applyFnc_ = applyObserver<o_t>;
  }

  // *** set or query tolerances ***

  // this is used to set a single tol for all
  template<class ToleranceType>
  void setTolerance(const ToleranceType & tolerance){
    correctionTolerances_.fill(tolerance);
    residualTolerances_.fill(tolerance);
    gradientTolerances_.fill(tolerance);
  }

  // finer-grained methods for tolerances
  void setCorrectionAbsoluteTolerance(correction_norm_type value){
    correctionTolerances_[0] = std::move(value); }
  void setCorrectionRelativeTolerance(correction_norm_type value){
    correctionTolerances_[1] = std::move(value); }

  void setResidualAbsoluteTolerance(residual_norm_type value) {
    residualTolerances_[0] = std::move(value); }
  void setResidualRelativeTolerance(residual_norm_type value) {
    residualTolerances_[1] = std::move(value); }

  void setGradientAbsoluteTolerance(gradient_norm_type value) {
    gradientTolerances_[0] = std::move(value); }
  void setGradientRelativeTolerance(gradient_norm_type value) {
    gradientTolerances_[1] = std::move(value); }

  auto correctionAbsoluteTolerance()const { return correctionTolerances_[0]; }
  auto correctionRelativeTolerance()const { return correctionTolerances_[1]; }
  auto residualAbsoluteTolerance()const   { return residualTolerances_[0]; }
  auto residualRelativeTolerance()const   { return residualTolerances_[1]; }
  auto gradientAbsoluteTolerance()const   { return gradientTolerances_[0]; }
  auto gradientRelativeTolerance()const   { return gradientTolerances_[1]; }

  template<typename SystemType>
  void solve(const SystemType & system, state_type & state)
  {
    // before we solve, we check if we need to recreate the updater
    // for example, this is the case if the system object changes
    if (recreateUpdater(system)){
      PRESSIOLOG_INFO("nonlinsolver: create updater");
      updater_ = createUpdater<this_type, SystemType>(state, updatingE_);
    }
    this->solveImpl(system, state, *updater_);
  }

#if not defined PRESSIO_ENABLE_TPL_PYBIND11
  template<class SystemType, class custom_updater_t>
  void solve(const SystemType & system,
	     state_type & state,
	     custom_updater_t && updater)
  {
    // here we have to create the updater each time becuase
    // we don't know if the updater changes or not

    PRESSIOLOG_INFO("nonlinsolver: custom updater");
    updatingE_ = ::pressio::nonlinearsolvers::Update::Custom;

    using u_t =
      mpl::conditional_t<
	std::is_lvalue_reference<custom_updater_t&&>::value,
      Updater<SystemType, state_type, this_type, custom_updater_t&>,
      Updater<SystemType, state_type, this_type, custom_updater_t>
      >;
    updater_ = pressio::utils::make_unique<u_t>(std::forward<custom_updater_t>(updater));
    updater_->applyFnc_ = applyUpdater<u_t>;
    updater_->resetFnc_ = resetUpdater<u_t>;
    this->solveImpl(system, state, *updater_);
  }

#else
  template<typename SystemType>
  void solveForPy(pybind11::object pySystem, state_type state)
  {
    SystemType system(pySystem);

    // before we solve, we check if we need to recreate the updater
    // for example, this is the case if the system object changes
    if (recreateUpdater(system)){
      PRESSIOLOG_INFO("nonlinsolver: create updater");
      updater_ = createUpdater<this_type, SystemType>(state, updatingE_);
    }
    this->solveImpl(system, state, *updater_);
  }
#endif

private:
  template<typename SystemType>
  bool recreateUpdater(const SystemType & system)
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

  template<class SystemType, class StateType, class UpdaterType>
  void solveImpl(const SystemType & system,
		 StateType & state,
		 UpdaterType & updater)
  {
    PRESSIOLOG_INFO("nonlinsolver: solve");

    residual_norm_type   residualNorm0 = {};
    correction_norm_type correctionNorm0 = {};
    gradient_norm_type   gradientNorm0 = {};
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
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e)
      {
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::NonlinearSolveFailure();
      }
      catch (::pressio::eh::ResidualHasNans const &e)
      {
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::NonlinearSolveFailure();
      }

      // 2.
      const auto correctionNorm = T::correctionNormCurrentCorrectionStep();
      const auto residualNorm	= T::residualNormCurrentCorrectionStep();
      if (iStep_==1) {
	residualNorm0   = residualNorm;
	correctionNorm0 = correctionNorm;
      }
      correctionNorms_[0] = correctionNorm;
      correctionNorms_[1] = correctionNorms_[0]/correctionNorm0;
      residualNorms_[0]   = residualNorm;
      residualNorms_[1]   = residualNorms_[0]/residualNorm0;

      if (T::hasGradientComputation()){
	const auto gradientNorm	= T::gradientNormCurrentCorrectionStep();
	if (iStep_==1) gradientNorm0 = gradientNorm;

	gradientNorms_[0] = gradientNorm;
	gradientNorms_[1] = gradientNorms_[0]/gradientNorm0;
      }

      if (T::hasGradientComputation()){
	impl::print_metrics
	  (iStep_, printStrippedMetrics_,
	   correctionNorms_[0], correctionNorms_[1],
	   residualNorms_[0],   residualNorms_[1],
	   gradientNorms_[0],   gradientNorms_[1]);
      }
      else{
	impl::print_metrics
	  (iStep_, printStrippedMetrics_,
	   correctionNorms_[0], correctionNorms_[1],
	   residualNorms_[0],   residualNorms_[1]);
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
  bool stopLoop(const iteration_type & iStep) const
  {
    switch (stoppingE_)
      {
      case Stop::AfterMaxIters:
    	return iStep == iterative_base_type::maxIters_;

      case Stop::WhenCorrectionAbsoluteNormBelowTolerance:
    	return correctionNorms_[0] < correctionTolerances_[0];
      case Stop::WhenCorrectionRelativeNormBelowTolerance:
    	return correctionNorms_[1] < correctionTolerances_[1];

      case Stop::WhenResidualAbsoluteNormBelowTolerance:
    	return residualNorms_[0] < residualTolerances_[0];
      case Stop::WhenResidualRelativeNormBelowTolerance:
    	return residualNorms_[1] < residualTolerances_[1];

      case Stop::WhenGradientAbsoluteNormBelowTolerance:
    	return gradientNorms_[0] < gradientTolerances_[0];
      case Stop::WhenGradientRelativeNormBelowTolerance:
    	return gradientNorms_[1] < gradientTolerances_[1];

      default:
    	return false;
      };
  }

  void fillTolerancesWithDefaultValues()
  {
    correctionTolerances_.fill(0.000001);
    residualTolerances_.fill(0.000001);
    gradientTolerances_.fill(0.000001);
  }

};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
