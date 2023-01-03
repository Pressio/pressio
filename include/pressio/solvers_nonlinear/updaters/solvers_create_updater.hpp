/*
//@HEADER
// ************************************************************************
//
// solvers_create_updater.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_CREATE_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_CREATE_UPDATER_HPP_

#include "solvers_updater.hpp"
#include "solvers_default.hpp"
// #include "solvers_armijo.hpp"
// #include "solvers_backtrack_strictly_decreasing_objective.hpp"
// #include "solvers_lm_schedule1.hpp"
// #include "solvers_lm_schedule2.hpp"

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template<typename T>
void apply_updater(BaseUpdater* updater,
		  const void* system_as_void,
		  void * state_as_void,
		  void * solver_as_void)
{
  using system_t = typename T::system_type;
  using state_t = typename T::state_type;
  using solver_t = typename T::solver_type;

  auto* p = static_cast<T*>(updater);
  auto* sys = reinterpret_cast<const system_t*>(system_as_void);
  auto* state = reinterpret_cast<state_t*>(state_as_void);
  auto* solver = reinterpret_cast<solver_t*>(solver_as_void);
  p->get()(*sys, *state, *solver);
}

template<typename T>
void reset_updater(BaseUpdater* updater)
{
  auto* p = static_cast<T*>(updater);
  p->get().reset();
}

/*
  GaussNewton only admits: standard, Armijo
*/
template<class SolverType, class SystemType, typename StateType>
mpl::enable_if_t<
  std::is_same<typename SolverType::tag_type, GaussNewton>::value or
  std::is_same<typename SolverType::tag_type, GaussNewtonQR>::value,
  std::unique_ptr<BaseUpdater>
  >
create_updater(const StateType & state,
	       ::pressio::nonlinearsolvers::Update updateE)
{
  using res_t = std::unique_ptr<BaseUpdater>;

  switch (updateE)
    {
    case Update::Standard:{
      using f_t = DefaultUpdater;
      using u_t = Updater<SystemType, StateType, SolverType, f_t>;
      res_t result = pressio::utils::make_unique<u_t>(f_t{});
      result->applyFnc_ = apply_updater<u_t>;
      result->resetFnc_ = reset_updater<u_t>;
      return result;
    }

    // case Update::Armijo:{
    //   using f_t = ArmijoUpdater<StateType>;
    //   f_t F(state);
    //   using u_t = Updater<SystemType, StateType, SolverType, f_t>;
    //   res_t result = pressio::utils::make_unique<u_t>(std::move(F));
    //   result->applyFnc_ = apply_updater<u_t>;
    //   result->resetFnc_ = reset_updater<u_t>;
    //   return result;
    // }

    // case Update::BacktrackStrictlyDecreasingObjective:{
    //   using f_t = BacktrackStrictlyDecreasingObjectiveUpdater<StateType>;
    //   f_t F(state);
    //   using u_t = Updater<SystemType, StateType, SolverType, f_t>;
    //   res_t result = pressio::utils::make_unique<u_t>(std::move(F));
    //   result->applyFnc_ = apply_updater<u_t>;
    //   result->resetFnc_ = reset_updater<u_t>;
    //   return result;
    // }

    default:
      throw std::runtime_error("Invalid update enum for GaussNewton");
      return nullptr;
    }
}

// /*
//    Levenberg-Mqrquardt: standard, LM1, LM2
// */
// template<typename solver_t, typename SystemType, typename StateType>
// mpl::enable_if_t<
//   std::is_same<typename solver_t::solver_tag, LevenbergMarquardt>::value,
//   std::unique_ptr<impl::BaseUpdater>
//   >
// createUpdater(const StateType & state,
// 	      ::pressio::nonlinearsolvers::Update updateE)
// {
//   using res_t = std::unique_ptr<BaseUpdater>;

//   switch (updateE)
//     {
//     case Update::Standard:{
//       using f_t = DefaultUpdater;
//       using u_t = Updater<SystemType, StateType, solver_t, f_t>;
//       res_t result = pressio::utils::make_unique<u_t>(f_t{});
//       result->applyFnc_ = apply_updater<u_t>;
//       result->resetFnc_ = reset_updater<u_t>;
//       return result;
//     }

//     case Update::LMSchedule1:{
//       using f_t = LMSchedule1Updater<StateType>;
//       f_t F(state);
//       using u_t = Updater<SystemType, StateType, solver_t, f_t>;
//       res_t result = pressio::utils::make_unique<u_t>(std::move(F));
//       result->applyFnc_ = apply_updater<u_t>;
//       result->resetFnc_ = reset_updater<u_t>;
//       return result;
//     }

//     case Update::LMSchedule2:{
//       using f_t = LMSchedule2Updater<StateType>;
//       f_t F(state);
//       using u_t = Updater<SystemType, StateType, solver_t, f_t>;
//       res_t result = pressio::utils::make_unique<u_t>(std::move(F));
//       result->applyFnc_ = apply_updater<u_t>;
//       result->resetFnc_ = reset_updater<u_t>;
//       return result;
//     }

//     default:
//       throw std::runtime_error("Invalid update enum for LevenbergMarquardt");
//       return nullptr;
//     }
// }

/*
   Newton-Raphson
*/
template<typename solver_t, typename SystemType, typename StateType>
mpl::enable_if_t<
  std::is_same<typename solver_t::tag_type, NewtonRaphson>::value,
  std::unique_ptr<impl::BaseUpdater>
  >
create_updater(const StateType & state,
	       ::pressio::nonlinearsolvers::Update updateE)
{
  using res_t = std::unique_ptr<BaseUpdater>;

  switch (updateE)
    {
    case Update::Standard:{
      using f_t = DefaultUpdater;
      using u_t = Updater<SystemType, StateType, solver_t, f_t>;
      res_t result = pressio::utils::make_unique<u_t>(f_t{});
      result->applyFnc_ = apply_updater<u_t>;
      result->resetFnc_ = reset_updater<u_t>;
      return result;
    }

    // case Update::BacktrackStrictlyDecreasingObjective:{
    //   using f_t = BacktrackStrictlyDecreasingObjectiveUpdater<StateType>;
    //   f_t F(state);
    //   using u_t = Updater<SystemType, StateType, solver_t, f_t>;
    //   res_t result = pressio::utils::make_unique<u_t>(std::move(F));
    //   result->applyFnc_ = apply_updater<u_t>;
    //   result->resetFnc_ = reset_updater<u_t>;
    //   return result;
    // }

    default:
      throw std::runtime_error("Invalid update enum for NewtonRaphson");
      return nullptr;
    }
}

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_CREATE_UPDATER_HPP_
