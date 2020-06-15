/*
//@HEADER
// ************************************************************************
//
// solvers_lm_neq_scheulde_policies.hpp
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

#ifndef SOLVERS_LM_NEQ_SCHEDULE_POLICIES_HPP_
#define SOLVERS_LM_NEQ_SCHEDULE_POLICIES_HPP_

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template<typename lm_schedule_policy_tag, typename scalar_t>
class LMSchedule;

/*
Updating technique from 
DAMPING PARAMETER IN MARQUARDT’S METHOD
Hans Bruun Nielsen
TECHNICAL REPORT IMM-REP-1999-05
Update given in Equation 2.5

Note that we solve J^T J + diag(J^TJ)\mu  here
*/
template<typename scalar_t>
class LMSchedule<pressio::solvers::iterative::lm::SchedulePolicyDefault, scalar_t>
{
private:
    const scalar_t beta_lm_;
    const scalar_t gamma_lm_;
    const scalar_t p_lm_;
    const scalar_t tau_lm_;
    scalar_t nu_lm_;
    scalar_t mu_ = 0.;

public:
    LMSchedule():
      beta_lm_{2.},
      gamma_lm_{3.},
      p_lm_{3.},
      nu_lm_{2.},
      tau_lm_{1.}
    {}
    LMSchedule(scalar_t beta_lm, 
               scalar_t gamma_lm, 
               scalar_t p_lm, 
               scalar_t tau_lm) : 
               beta_lm_(beta_lm),
               gamma_lm_(gamma_lm),
               p_lm_(p_lm),
               nu_lm_(beta_lm),
               tau_lm_(tau_lm)
    {}
   
    scalar_t getMu(){
      return mu_;
    }

    template<typename hessian_t>
    void reset(hessian_t H){
      mu_ = tau_lm_;
    }

    
    template <
      typename system_t,
      typename gradient_t,
      typename residual_t,
      typename hessian_t,
      typename state_t
      >
    void evaluate(state_t & stateInOut, 
                  state_t & ytrial,
                  const gradient_t & correction,
                  const gradient_t & gradient, 
                  residual_t & residual, 
                  const hessian_t & hessian, 
                  const system_t & sys)
    {
      constexpr auto one = ::pressio::utils::constants<scalar_t>::one();
      ::pressio::ops::do_update(ytrial, stateInOut, one, correction,one);
      auto & HObj = *hessian.data();
      gradient_t tmpa(correction.extent(0));
      for (int i=0; i< hessian.extent(0); i++){
        (*tmpa.data())(i) = (*correction.data())(i)*HObj(i,i);
      } 
      auto hh = ::pressio::ops::dot(correction,tmpa);
      auto hg = ::pressio::ops::dot(correction,gradient);
      auto denom = 0.5*(mu_*hh + hg);  //note sign difference in gradient
      auto r2_old = ::pressio::ops::norm2(residual); 
      r2_old = r2_old*r2_old;
      sys.residual(ytrial, residual);
      auto r2_new = ::pressio::ops::norm2(residual);
      r2_new = r2_new*r2_new; 
      auto rho = (r2_old - r2_new) / (denom);
      /*
      This block of code is used if we solve J^T J + I mu = J^TR
      auto hh = ::pressio::ops::dot(correction,correction);
      auto hg = ::pressio::ops::dot(correction,gradient);
      auto denom = 0.5*(mu_*hh + hg);  //note sign difference in gradient
      auto r2_old = ::pressio::ops::dot(residual,residual); 
      sys.residual(ytrial, residual);
      auto r2_new = ::pressio::ops::dot(residual,residual); 
      auto rho = (r2_old - r2_new) / (denom);
      */
      if (rho > 0){
        ::pressio::ops::do_update(stateInOut, one, correction, one);
        scalar_t mu_rat = 1 - (beta_lm_ - 1.)*pow(2.*rho-1,p_lm_);
        mu_rat = std::max(mu_rat,1./gamma_lm_);
        mu_ = mu_*mu_rat;
        nu_lm_ = beta_lm_;
      } else{
        mu_ = mu_*nu_lm_;
        nu_lm_ = 2.*nu_lm_; 
      }
//      ::pressio::utils::io::print_stdout(std::scientific,
//  				    "mu =", mu_,
//  				    utils::io::reset(),
//  				    "\n");
    } 

};

/*
Updating technique from 
DAMPING PARAMETER IN MARQUARDT’S METHOD
Hans Bruun Nielsen
TECHNICAL REPORT IMM-REP-1999-05

Update given in Equation 2.4
*/
template<typename scalar_t>
class LMSchedule<pressio::solvers::iterative::lm::SchedulePolicy2,scalar_t>
{

private:
    const scalar_t rho1_lm_;
    const scalar_t rho2_lm_;
    const scalar_t beta_lm_;
    const scalar_t gamma_lm_;
    const scalar_t tau_lm_;
    scalar_t mu_ = 0.;
public:
    LMSchedule():
      rho1_lm_{0.2},
      rho2_lm_{0.8},
      beta_lm_{2.},
      gamma_lm_{3.},
      tau_lm_{1.}
    {}
    LMSchedule(scalar_t beta_lm,
               scalar_t gamma_lm, 
               scalar_t rho1_lm, 
               scalar_t rho2_lm, 
               scalar_t tau_lm) : 
               beta_lm_(beta_lm),
               gamma_lm_(gamma_lm),
               rho1_lm_(rho1_lm),
               rho2_lm_(rho2_lm),
               tau_lm_(tau_lm)
    {}

    scalar_t getMu(){
      return mu_;
    }

    template<typename hessian_t>
    void reset(hessian_t H){
      mu_ = tau_lm_;
    }


    template <
      typename system_t,
      typename gradient_t,
      typename residual_t,
      typename hessian_t,
      typename state_t
      >
    void evaluate(state_t & stateInOut, 
                  state_t & ytrial,
                  const gradient_t & correction,
                  const gradient_t & gradient, 
                  residual_t & residual, 
                  const hessian_t & hessian, 
                  const system_t & sys)
    {
      constexpr auto one = ::pressio::utils::constants<scalar_t>::one();
      ::pressio::ops::do_update(ytrial, stateInOut, one, correction,one);
      auto & HObj = *hessian.data();
      gradient_t tmpa(correction.extent(0));
      for (int i=0; i< hessian.extent(0); i++){
        (*tmpa.data())(i) = (*correction.data())(i)*HObj(i,i);
      } 
      auto hh = ::pressio::ops::dot(correction,tmpa);
      auto hg = ::pressio::ops::dot(correction,gradient);
      auto denom = 0.5*(mu_*hh + hg);  //note sign difference in gradient
      auto r2_old = ::pressio::ops::norm2(residual);
      r2_old = r2_old*r2_old; 
      sys.residual(ytrial, residual);
      auto r2_new = ::pressio::ops::norm2(residual);
      r2_new = r2_new*r2_new; 
      auto rho = (r2_old - r2_new) / (denom);
      if (rho < rho1_lm_){
        mu_ = std::min(mu_*beta_lm_,pow(10.,7.));
      } 
      if (rho > rho2_lm_){
        mu_ = std::max(pow(10.,-7.),mu_/gamma_lm_);
      }
      if (rho > 0){
        ::pressio::ops::do_update(stateInOut, one, correction,one);
      }
//      ::pressio::utils::io::print_stdout(std::scientific,
//  				    "mu =", mu_,
//  				    utils::io::reset(),
//  				    "\n");
    } 

};



}}}} //end namespace pressio::solvers::iterative::impl
#endif
