/*
//@HEADER
// ************************************************************************
//
// rom_solve_problem_functions.hpp
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

#ifndef ROM_LSPG_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_
#define ROM_LSPG_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_

#include "./impl/unsteady/rom_unsteady_problem_solve_functions_impl.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/* these functions make it easier for the users
   so that they pass a rom problem and don't need to
   worry about knowing that they need to get the stepper
   and pass that to the ode integrators or the system to
   the solver for the steady case */

/************************************
	      STEADY
***********************************/
template<typename rom_problem_t, typename rom_state_t, typename solver_t>
void solveSteady(rom_problem_t & problem,
		 rom_state_t & romState,
		 solver_t & solver)
{
  static_assert
    (::pressio::rom::details::traits<rom_problem_t>::is_steady_lspg,
     "rom::lspg::solve() can only be called for a steady problem");

  solver.solve(problem.systemRef(), romState);
}

/************************************
	      UNSTEADY
***********************************/
/*--------------------------------------------
  solve for fixed number of steps
  --------------------------------------------*/
template<typename rom_problem_type, typename ...Args>
void solveNTimes
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}

// For pressio4py, I cannot handle the variadic case so I need to specify
// the number of arguments and this will work. For now leave to just
// a specific case, we can add more later.
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
// with valid collector
template<
  class rom_problem_type,
  class rom_state_t,
  class timet,
  class collector_t,
  class solver_t
  >
void solveNSequentialMinimizations(rom_problem_type & problem,
				   typename rom_state_t::traits::wrapped_t & stateInOut,
				   timet t0,
				   timet dt,
				   const ::pressio::ode::types::step_t numSteps,
				   pybind11::object pyCollector,
				   solver_t & solver)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  rom_state_t stateView(stateInOut, ::pressio::view());

  collector_t collector(pyCollector);
  ::pressio::ode::advanceNSteps
      (problem.stepperRef(), stateView, t0, dt, numSteps, collector, solver);
}

// without collector
template<
  class rom_problem_type, class rom_state_t, class timet, class solver_t
  >
void solveNSequentialMinimizations(rom_problem_type & problem,
				   typename rom_state_t::traits::wrapped_t & stateInOut,
				   timet t0,
				   timet dt,
				   const ::pressio::ode::types::step_t numSteps,
				   solver_t & solver)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  rom_state_t stateView(stateInOut, ::pressio::view());

  ::pressio::ode::advanceNSteps
      (problem.stepperRef(), stateView, t0, dt, numSteps, solver);
}
#else

template<typename rom_problem_type, typename ...Args>
void solveNSequentialMinimizations
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}
#endif


template<typename rom_problem_type, typename ...Args>
void solveNSequentialResidualMinimizations
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialResidualMinimizationProblemNTimes
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}

/*--------------------------------------------
  solve to target time
  -------------------------------------------- */
template<typename rom_problem_type, typename ...Args>
void solveToTargetTime
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTime(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialMinimizationsToTargetTime
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTime(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialResidualMinimizationsToTargetTime
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTime(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialResidualMinimizationProblemToTargetTime
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTime(problem, std::forward<Args>(args)...);
}

/*--------------------------------------------
  advance to target time with step recovery
  -------------------------------------------- */
template<typename rom_problem_type, typename ...Args>
void solveToTargetTimeWithRecovery
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTimeWithRec(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialMinimizationsToTargetTimeWithRecovery
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTimeWithRec(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialResidualMinimizationsToTargetTimeWithRecovery
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTimeWithRec(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveSequentialResidualMinimizationProblemToTargetTimeWithRecovery
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyToTimeWithRec(problem, std::forward<Args>(args)...);
}

}}}//end namespace pressio::rom::lspg
#endif  // ROM_LSPG_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_
