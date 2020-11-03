/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_solve_functions.hpp
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

#ifndef ROM_LSPG_UNSTEADY_PROBLEM_SOLVE_FUNCTIONS_HPP_
#define ROM_LSPG_UNSTEADY_PROBLEM_SOLVE_FUNCTIONS_HPP_

#include "./impl/rom_lspg_unsteady_problem_solve_functions_impl.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/* these are here just to make it easier for the users
   so that they pass a rom problem and don't need to
   worry about knowing that they need to get the stepper
   and pass that to the ode integrators or the system to
   the solver for the steady case */

/************************************
	      STEADY
***********************************/
// allow romState and solver to be passed in both orders
template<typename rom_problem_t, typename rom_state_t, typename solver_t>
mpl::enable_if_t<
  ::pressio::rom::concepts::rom_state<rom_state_t>::value and
  !::pressio::rom::concepts::rom_state<solver_t>::value
  >
solveSteady(rom_problem_t & problem, rom_state_t & romState, solver_t & solver)
{
  static_assert
    (::pressio::rom::details::traits<rom_problem_t>::is_steady_lspg,
     "rom::lspg::solve() can only be called for a steady problem");

  solver.solve(problem.systemRef(), romState);
}

template<typename rom_problem_t, typename rom_state_t, typename solver_t>
mpl::enable_if_t<
  ::pressio::rom::concepts::rom_state<rom_state_t>::value and
  !::pressio::rom::concepts::rom_state<solver_t>::value
  >
solveSteady(rom_problem_t & problem, solver_t & solver, rom_state_t & romState)
{
  solveSteady(problem, romState, solver);
}


/************************************
	      UNSTEADY
***********************************/

/* fixed number of steps */
template<typename rom_problem_type, typename ...Args>
void solveNTimes
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}

template<typename rom_problem_type, typename ...Args>
void solveNSequentialMinimizations
(rom_problem_type & problem, Args && ...args)
{
  impl::_lspgUnsteadyNTimes(problem, std::forward<Args>(args)...);
}

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

/* to target time */
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

/* to target time with step recovery */
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
#endif  // ROM_ROM_UNSTEADY_PROBLEM_ADVANCERS_HPP_
