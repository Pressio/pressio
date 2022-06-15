/*
//@HEADER
// ************************************************************************
//
// rom_entropy.hpp
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

#ifndef PRESSIO_ROM_GALERKIN_ENTROPY_TOPLEVEL_INC_HPP_
#define PRESSIO_ROM_GALERKIN_ENTROPY_TOPLEVEL_INC_HPP_

#include "./mpl.hpp"
#include "./utils.hpp"
#include "./type_traits.hpp"
#include "./ops.hpp"
#include "./qr.hpp"

#include "./rom/predicates/all.hpp"
#include "./rom/constraints/all.hpp"
#include "./rom_decoder.hpp"
#include "./rom/rom_reconstructor_fom_state.hpp"
#include "./rom/rom_manager_fom_states.hpp"
//
//
//
//
namespace pressio{ namespace rom{ namespace galerkin_entropy{

template <typename fom_t, typename rom_state_t,  typename basis_t> 
class galerkinEntropyRom
{
  rom_state_t & romState_;
  const fom_t & appObj_;
  const basis_t & Phi_;
public:
  galerkinEntropyRom(const fom_t & appObj, rom_state_t & romState,const basis_t & Phi) : appObj_(appObj), romState_(romState), Phi_(Phi)  {};

  template <typename obs_t, typename scalar_t>
  void advance_n_steps_rk4(rom_state_t romState,int nSteps, scalar_t dt,scalar_t tau, obs_t & Obs){
    auto reducedVelocity = ::pressio::ops::clone(romState);
    auto fomState = appObj_.createVelocity();
    auto fomVelocity = appObj_.createVelocity();
    auto fomProjectedVelocity = appObj_.createVelocity();
    auto MPhi = ::pressio::ops::clone(Phi_);
    auto K = romState.size();
    Eigen::MatrixXd reducedMassEigen(K,K);
    Kokkos::View<scalar_t**,Kokkos::LayoutLeft> reducedMass("reducedMass",K,K);
    rom_state_t MinvF(K);

    scalar_t rk4const[4];
    rk4const[0] = 1./4.;
    rk4const[1] = 1./3.;
    rk4const[2] = 1./2.;
    rk4const[3] = 1.;
    
    scalar_t epsilon = 1.e-4; 
    for (int step = 0; step <= nSteps; step++){
      auto romState0 = ::pressio::ops::clone(romState);
      for (int k = 0; k < 4 ; k++){
        //reconstruct state
        ::pressio::ops::product(::pressio::nontranspose(), 1., Phi_,romState,0.,fomState);
        // evaluate the velocity
        auto timeOne = std::chrono::high_resolution_clock::now();
        scalar_t timeAtStep = dt*step;
        appObj_.velocity(fomState,timeAtStep,fomVelocity);

        auto timeTwo = std::chrono::high_resolution_clock::now(); 
        const std::chrono::duration<double> elapsedOne = timeTwo - timeOne;
       // project the velocity
       ::pressio::ops::product(::pressio::transpose(),1., Phi_,fomVelocity, 0.,reducedVelocity);    

       auto timeThree = std::chrono::high_resolution_clock::now(); 
       const std::chrono::duration<double> elapsedTwo = timeThree - timeTwo;

       // compute the mass matrix and project
       appObj_.applyMassMatrix(fomState,Phi_,MPhi);

       auto timeFour = std::chrono::high_resolution_clock::now(); 
       const std::chrono::duration<double> elapsedThree = timeFour - timeThree;

       ::pressio::ops::product(::pressio::transpose(),::pressio::nontranspose(),1., Phi_,MPhi, 0.,reducedMass);    
       for (int n=0; n < K; n++){
         for (int m=0; m< K ; m++){
           reducedMassEigen(m,n) = reducedMass(m,n);
         }
       } 

       auto timeFive = std::chrono::high_resolution_clock::now(); 
       const std::chrono::duration<double> elapsedFour = timeFive - timeFour;

       // Now we need to solve the linear system  for the velocity
       //MinvF = reducedMassEigen.colPivHouseholderQr().solve(reducedVelocity);
       MinvF = reducedMassEigen.householderQr().solve(reducedVelocity);

       auto timeSix = std::chrono::high_resolution_clock::now(); 
       const std::chrono::duration<double> elapsedFive = timeSix - timeFive;

       // step in time
       ::pressio::ops::update(romState,0.,romState0,1.,MinvF,dt*rk4const[k]);
      std::cout << " velocity time = " << elapsedOne.count() << std::endl;
      std::cout << " V^T f time = " << elapsedTwo.count() << std::endl;
      std::cout << " MPhi time = " << elapsedThree.count() << std::endl;
      std::cout << " Phi x MPhi time = " << elapsedFour.count() << std::endl;
      std::cout << " solve time = " << elapsedFive.count() << std::endl;

      }
      Obs(step,dt*step,romState);
      std::cout << " Step = " << step << std::endl;

    }
  }

};

}}} // end name space 

#endif
