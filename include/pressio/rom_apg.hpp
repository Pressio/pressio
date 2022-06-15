/*
//@HEADER
// ************************************************************************
//
// rom_apg.hpp
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

#ifndef PRESSIO_ROM_APG_TOPLEVEL_INC_HPP_
#define PRESSIO_ROM_APG_TOPLEVEL_INC_HPP_

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
//#include "./rom/rom_public_api_apg.hpp"
//
//
//
//
namespace pressio{ namespace rom{ namespace apg{

template <typename fom_t, typename rom_state_t,  typename basis_t> 
class apgRom
{
  rom_state_t & romState_;
  const fom_t & appObj_;
  const basis_t & Phi_;
public:
  apgRom(const fom_t & appObj, rom_state_t & romState,const basis_t & Phi) : appObj_(appObj), romState_(romState), Phi_(Phi)  {};



  template <typename obs_t, typename scalar_t>
  void advance_n_steps_rk4_with_mass(rom_state_t & romState,int nSteps, scalar_t dt,scalar_t tau, obs_t & Obs){
    auto reducedVelocity = ::pressio::ops::clone(romState);
    auto fomState = appObj_.createVelocity();
    auto fomVelocity = appObj_.createVelocity();
    auto MInvTimesFomVelocity = appObj_.createVelocity();
    auto fomProjectedVelocity = appObj_.createVelocity();
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
        scalar_t timeAtStep = dt*step;
        appObj_.velocity(fomState,timeAtStep,fomVelocity);
        
       // project the velocity
       ::pressio::ops::product(::pressio::transpose(),1., Phi_,fomVelocity, 0.,reducedVelocity);    

       // compute the orthogonal compliment to the   velocity,(M^{-1} -  \Phi \Phi^T) f
       //                                                    = M^{-1}f - \Phi\Phi^T f
       ::pressio::ops::product(::pressio::nontranspose(), 1., Phi_,reducedVelocity,0.,fomProjectedVelocity);
       ::pressio::ops::update(MInvTimesFomVelocity,0.,fomVelocity,1.);
       appObj_.applyInverseMassMatrix(MInvTimesFomVelocity);
       ::pressio::ops::update(fomProjectedVelocity,0.,MInvTimesFomVelocity,1.,fomProjectedVelocity,-1.);

       // Compute matrix free approximation to J (I - \Pi) f
       // firs tcompute f( x + eps (I - \PI f) )
       ::pressio::ops::update(fomState,0.,fomState,1.,fomProjectedVelocity,epsilon);
       appObj_.velocity(fomState,timeAtStep,fomProjectedVelocity);

       // matrix free approximation, Jf = f(x + eps f) - f(x) / eps                
       ::pressio::ops::update(fomProjectedVelocity,0.,fomProjectedVelocity,1./epsilon,fomVelocity,-1./epsilon);

       // project the stabilization term  and add to reduced velocity
       ::pressio::ops::product(::pressio::transpose(),tau, Phi_,fomProjectedVelocity, 1.,reducedVelocity);    

       // step in time
       ::pressio::ops::update(romState,0.,romState0,1.,reducedVelocity,dt*rk4const[k]);
      }
      Obs(step,dt*step,romState);
//      std::cout << " Step = " << step << std::endl;
//      std::cout << " Norm = " << reducedVelocity.norm() << std::endl;
    }
  }




  template <typename obs_t, typename scalar_t>
  void advance_n_steps_rk4(rom_state_t romState,int nSteps, scalar_t dt,scalar_t tau, obs_t & Obs){
    auto reducedVelocity = ::pressio::ops::clone(romState);
    auto fomState = appObj_.createVelocity();
    auto fomVelocity = appObj_.createVelocity();
    auto fomProjectedVelocity = appObj_.createVelocity();
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
        scalar_t timeAtStep = dt*step;
        appObj_.velocity(fomState,timeAtStep,fomVelocity);
        
       // project the velocity
       ::pressio::ops::product(::pressio::transpose(),1., Phi_,fomVelocity, 0.,reducedVelocity);    

       // compute the orthogonal complliment to the   velocity,(I -  \Phi \Phi^T) f
       ::pressio::ops::product(::pressio::nontranspose(), 1., Phi_,reducedVelocity,0.,fomProjectedVelocity);
       ::pressio::ops::update(fomProjectedVelocity,0.,fomVelocity ,1.,fomProjectedVelocity,-1.);

       // Compute matrix free approximation to J (I - \Pi) f
       // firs tcompute f( x + eps (I - \PI f) )
       ::pressio::ops::update(fomState,0.,fomState,1.,fomProjectedVelocity,epsilon);
       appObj_.velocity(fomState,timeAtStep,fomProjectedVelocity);

       // matrix free approximation, Jf = f(x + eps f) - f(x) / eps                
       ::pressio::ops::update(fomProjectedVelocity,0.,fomProjectedVelocity,1./epsilon,fomVelocity,-1./epsilon);

       // project the stabilization term  and add to reduced velocity
       ::pressio::ops::product(::pressio::transpose(),tau, Phi_,fomProjectedVelocity, 1.,reducedVelocity);    

       // step in time
       ::pressio::ops::update(romState,0.,romState0,1.,reducedVelocity,dt*rk4const[k]);
      }
      Obs(step,dt*step,romState);
      std::cout << " Step = " << step << std::endl;
    }
  }

};

}}} // end name space apg

#endif
