/*
//@HEADER
// ************************************************************************
//
// solvers_apply_updater.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_APPLY_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_APPLY_UPDATER_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

/*
  GaussNewton only admits: standard, armijo
*/
template<
  typename system_t,
  typename state_t,
  typename solver_t
  >
mpl::enable_if_t<
  std::is_same<typename solver_t::solver_tag, GaussNewton>::value or
  std::is_same<typename solver_t::solver_tag, GaussNewtonQR>::value
  >
applyUpdater(const system_t & system,
	     state_t & state,
	     solver_t & solver,
	     ::pressio::solvers::nonlinear::update updateE,
	     std::shared_ptr<impl::BaseUpdater> updater)
{
  switch (updateE)
  {
    case update::standard:{
      using ut = impl::DefaultUpdater;
      updater->template updateState<ut>(system, state, solver);
      break;
    }

    case update::armijo:{
      using ut = impl::ArmijoUpdater<state_t>;
      updater->template updateState<ut>(system, state, solver);
      break;
    }

    default:
      throw std::runtime_error("Invalid update enum for GaussNewton");
  }
}


/*
  Levenberg-Mqrquardt: standard, LM1, LM2
*/
template<
  typename system_t,
  typename state_t,
  typename solver_t
  >
mpl::enable_if_t<
  std::is_same<typename solver_t::solver_tag, LevenbergMarquardt>::value
  >
applyUpdater(const system_t & system,
	     state_t & state,
	     solver_t & solver,
	     ::pressio::solvers::nonlinear::update updateE,
	     std::shared_ptr<impl::BaseUpdater> updater)
{
  switch (updateE)
  {
    case update::standard:{
      using ut = impl::DefaultUpdater;
      updater->template updateState<ut>(system, state, solver);
      break;
    }

    case update::LMSchedule1:{
	    using ut = impl::LMSchedule1Updater<state_t>;
	    updater->template updateState<ut>(system, state, solver);
	    break;
    }

    case update::LMSchedule2:{
	    using ut = impl::LMSchedule2Updater<state_t>;
	    updater->template updateState<ut>(system, state, solver);
	    break;
    }

    default:
      throw std::runtime_error("Invalid update enum for LevenbergMarquardt");
  }
}

/*
  Newton-Raphson admits: standard
*/
template<
  typename system_t,
  typename state_t,
  typename solver_t
  >
mpl::enable_if_t<
  std::is_same<typename solver_t::solver_tag, NewtonRaphson>::value
  >
applyUpdater(const system_t & system,
	     state_t & state,
	     solver_t & solver,
	     ::pressio::solvers::nonlinear::update updateE,
	     std::shared_ptr<impl::BaseUpdater> updater)
{
  switch (updateE)
  {
    case update::standard:{
      using ut = impl::DefaultUpdater;
      updater->template updateState<ut>(system, state, solver);
      break;
    }

    default:
      throw std::runtime_error("Invalid update enum for NewtonRaphson");
  }
}

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_APPLY_UPDATER_HPP_
