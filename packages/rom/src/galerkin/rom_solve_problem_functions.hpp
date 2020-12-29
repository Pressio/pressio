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

#ifndef ROM_GALERKIN_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_
#define ROM_GALERKIN_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_

namespace pressio{ namespace rom{ namespace galerkin{

/* these functions make it easier for the users because
   they pass a rom problem and don't need to extract from
   the problem the stepper and pass that to the integrators */

/*--------------------------------------------
  solve for fixed number of steps
  -------------------------------------------- */
// for pressio4py I cannot use variadic templates,
// so need to specify the various cases, for now
// just expose some subcases, not all

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
// explicit with collector
template<class rom_problem_type, class state_t, class time_t, class collector_t>
mpl::enable_if_t<rom_problem_type::stepper_t::is_explicit>
solveNSteps(rom_problem_type & problem,
	    typename state_t::traits::wrapped_t & stateInOut,
	    time_t t0, time_t dt,
	    ::pressio::ode::types::step_t num_steps,
	    pybind11::object pyCollector)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  state_t stateView(stateInOut, ::pressio::view());
  collector_t collector(pyCollector);
  ::pressio::ode::advanceNSteps
      (problem.stepperRef(), stateView, t0, dt, num_steps, collector);
}

// explicit wihtout collector
template<class rom_problem_type, class state_t, class time_t>
mpl::enable_if_t<rom_problem_type::stepper_t::is_explicit>
solveNSteps(rom_problem_type & problem,
	    typename state_t::traits::wrapped_t & stateInOut,
	    time_t t0, time_t dt,
	    ::pressio::ode::types::step_t num_steps)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  state_t stateView(stateInOut, ::pressio::view());
  ::pressio::ode::advanceNSteps
    (problem.stepperRef(), stateView, t0, dt, num_steps);
}

// implicit with collector
template<
  class rom_problem_type, class state_t, class time_t,
  class collector_t, class solver_t
  >
mpl::enable_if_t<rom_problem_type::stepper_t::is_implicit>
solveNSteps(rom_problem_type & problem,
	    typename state_t::traits::wrapped_t & stateInOut,
	    time_t t0, time_t dt,
	    ::pressio::ode::types::step_t num_steps,
	    pybind11::object pyCollector,
	    solver_t & solver)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  state_t stateView(stateInOut, ::pressio::view());
  collector_t collector(pyCollector);
  ::pressio::ode::advanceNSteps
      (problem.stepperRef(), stateView, t0, dt, num_steps, collector, solver);
}

// implicit without collector
template<class rom_problem_type, class state_t, class time_t, class solver_t>
mpl::enable_if_t<rom_problem_type::stepper_t::is_implicit>
solveNSteps(rom_problem_type & problem,
	    typename state_t::traits::wrapped_t & stateInOut,
	    time_t t0, time_t dt,
	    ::pressio::ode::types::step_t num_steps,
	    solver_t & solver)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  state_t stateView(stateInOut, ::pressio::view());
  ::pressio::ode::advanceNSteps
      (problem.stepperRef(), stateView, t0, dt, num_steps, solver);
}

#else

template<typename rom_problem_type, typename ...Args>
void solveNSteps
(rom_problem_type & problem, Args && ...args)
{
  ::pressio::ode::advanceNSteps
    (problem.stepperRef(), std::forward<Args>(args)...);
}
#endif

/*--------------------------------------------
  solve to target time
  -------------------------------------------- */
template<typename rom_problem_type, typename ...Args>
void solveToTargetTime
(rom_problem_type & problem, Args && ...args)
{
  ::pressio::ode::advanceToTargetTime
    (problem.stepperRef(), std::forward<Args>(args)...);
}

/*--------------------------------------------
  solve to target time with step recovery
  -------------------------------------------- */
template<typename rom_problem_type, typename ...Args>
void solveToTargetTimeWithTimeStepRecovery
(rom_problem_type & problem, Args && ...args)
{
  ::pressio::ode::advanceToTargetTimeWithTimeStepRecovery
    (problem.stepperRef(), std::forward<Args>(args)...);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_GALERKIN_ROM_SOLVE_PROBLEM_FUNCTIONS_HPP_
