/*
//@HEADER
// ************************************************************************
//
// rom_wls_hessian_gradient_system_api.hpp
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

#ifndef ROM_WLS_ROM_WLS_ADVANCE_ONE_WINDOW_HPP_
#define ROM_WLS_ROM_WLS_ADVANCE_ONE_WINDOW_HPP_

namespace pressio{ namespace rom{ namespace wls{

template <
  typename wls_system_t,
  typename wls_state_type,
  typename solverType,
  typename scalar_type
  >
void advanceOneWindow(wls_system_t & wlsSystem,
		      wls_state_type & wlsState,
		      solverType & solver,
		      const window_size_t & windowIndex,
		      scalar_type dt)
{
  const auto numStepsInWindow = wlsSystem.numStepsInWindow();
  const auto romSize = wlsSystem.romSize();
  const auto timeStencilSize = wlsSystem.timeStencilSize();

  //const auto activeWindowIndex_ = windowIndex;
  const auto windowStartTime_	= windowIndex*dt*numStepsInWindow;
  const auto step_s_		= windowIndex*numStepsInWindow;
  wlsSystem.setTimeStepSize(dt);
  wlsSystem.setWindowStartTime(windowStartTime_);
  wlsSystem.setStepS(step_s_);

  auto & wlsStateIC = wlsSystem.wlsStateInitConditionRef();

  // set initial guess over window (needed to avoid bad initial guesses that yield NaN)
  for (window_size_t i = 0; i < numStepsInWindow; i++){
    auto wlsViewAssign = ::pressio::containers::span(wlsState, i*romSize, romSize);
    auto wlsViewCopy = ::pressio::containers::span(wlsStateIC, (timeStencilSize - 1)*romSize, romSize);
    ::pressio::ops::deep_copy(wlsViewAssign, wlsViewCopy);
  }

  // solve system
  solver.solve(wlsSystem, wlsState);

  // Loop to update the the wlsStateIC vector.
  // If we add multistep explicit methods, need to add things here.
  const auto start = std::max(0, timeStencilSize - numStepsInWindow);
  for (window_size_t i = 0; i < start; ++i)
    {
      auto wlsTmpState	      = ::pressio::containers::span(wlsStateIC, i*romSize,     romSize);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC, (i+1)*romSize, romSize);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }

  for (window_size_t i = start ; i < timeStencilSize; ++i)
    {
      auto wlsTmpState  = ::pressio::containers::span(wlsStateIC, i*romSize, romSize);
      const auto wlsTmpState2 = ::pressio::containers::span(wlsState,
							    (numStepsInWindow-timeStencilSize+i)*romSize,
							    romSize);
      ::pressio::ops::deep_copy(wlsTmpState, wlsTmpState2);
    }

  PRESSIOLOG_INFO("Wls: completed window ", windowIndex);
}// end advanceOneWindow

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <
  typename wls_system_t,
  typename wls_state_type,
  typename solverType,
  typename scalar_type
  >
void solveWindowsSequentially(wls_system_t & wlsSystem,
			      typename wls_state_type::traits::wrapped_t & wlsStateIn,
			      solverType & solver,
			      const window_size_t & numWindows,
			      scalar_type dt)
{
  // here we want to view the state since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  wls_state_type stateView(wlsStateIn, ::pressio::view());

  for (auto iWind = 0; iWind < numWindows; iWind++){
    ::pressio::rom::wls::advanceOneWindow(wlsSystem, stateView, solver, iWind, dt);
  }
}

#else

template <
  typename wls_system_t,
  typename wls_state_type,
  typename solverType,
  typename scalar_type
  >
void solveWindowsSequentially(wls_system_t & wlsSystem,
			      wls_state_type & wlsState,
			      solverType & solver,
			      const window_size_t & numWindows,
			      scalar_type dt)
{
  for (auto iWind = 0; iWind < numWindows; iWind++){
    ::pressio::rom::wls::advanceOneWindow(wlsSystem, wlsState, solver, iWind, dt);
  }
}

#endif

}}}
#endif  // ROM_WLS_ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_HPP_
