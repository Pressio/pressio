/*
//@HEADER
// ************************************************************************
//
// rom_utils_set_gen_coordinates.hpp
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

#ifndef ROM_UTILS_ROM_UTILS_SET_GEN_COORDINATES_HPP_
#define ROM_UTILS_ROM_UTILS_SET_GEN_COORDINATES_HPP_

namespace pressio{ namespace rom{ namespace utils{

template <
  typename scalar_t,
  typename linear_solver_t,
  typename basis_t,
  typename fom_state_t,
  typename rom_state_t>
void set_gen_coordinates_L2_projection(linear_solver_t & linearSolver,
				       const basis_t & phi,
				       const fom_state_t & fomState,
				       const fom_state_t & fomNominalState,
				       rom_state_t & romState)
{
  /* Compute the ROM coefficients from optimal L^2 projection of yFOM */
  using hessian_t	= typename linear_solver_t::matrix_type;
  constexpr auto zero   = ::pressio::utils::constants<scalar_t>::zero();
  constexpr auto one    = ::pressio::utils::constants<scalar_t>::one();
  constexpr auto negOne = ::pressio::utils::constants<scalar_t>::negOne();

  const auto romSize = romState.extent(0);

  //compute hessian for phi^T phi
  hessian_t H(romSize,romSize);
  ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), one, phi, zero, H);

  //create a vector to store fomState - fomNominalState
  fom_state_t b(fomState);
  pressio::ops::update(b, one, fomNominalState, negOne);

  // compute phi^T b
  rom_state_t r(romSize);
  ::pressio::ops::product(::pressio::transpose(), one, phi, b, zero, r);

  // solve system for optimal L2 projection
  // allow the hessian to be overwritten since H is local so it does not matter
  linearSolver.solveAllowMatOverwrite(H, r, romState);
}

}}} //end namespace pressio::rom::utils
#endif  // ROM_UTILS_ROM_UTILS_SET_GEN_COORDINATES_HPP_
