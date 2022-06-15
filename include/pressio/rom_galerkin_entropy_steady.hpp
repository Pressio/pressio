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

#ifndef PRESSIO_ROM_GALERKIN_ENTROPY_STEADY_TOPLEVEL_INC_HPP_
#define PRESSIO_ROM_GALERKIN_ENTROPY_STEADY_TOPLEVEL_INC_HPP_

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




/*
template <>
void gmres(matrix_op_t Af, vec_t b, vec_t x0, scalar_t tol, int maxiter_outer, int maxiter)
{
    int k_outer = 0
    auto bnorm = ::pressio::ops::norm2(b);
    scalar_t error = 1.
    auto r = ::pressio::ops::clone(x0);
    while (k_outer < maxiter_outer && error >= tol)
    {
      Af(x0,r);
      ::pressio::ops::update(r,-1.,b,1.);
      //r = b - Af(x0)
      //if (main.mpi_rank == 0 and printnorm==1):
      //  print('Outer true norm = ' + str(np.linalg.norm(r)))
      //cs = np.zeros(maxiter) #should be the same on all procs
      //sn = np.zeros(maxiter) #same on all procs
      auto e1 = ::pressio::ops::clone(b);
      ::pressio::ops::set_zero(e1);
      e1[0] = 1;

      auto rnorm = ::pressio::ops::norm2(r);
      Q = np.zeros((np.size(b),maxiter)) 

      auto Q0 = ::pressio::ops::clone(r);
      ::pressio::ops::scale(Q0,rnorm);
      Q[:,0] = r / rnorm ## The first index of Q is across all procs
      H = np.zeros((maxiter + 1, maxiter)) ### this should be the same on all procs
      auto beta = ::pressio::ops::clone(e1);
      ::pressio::ops::scale(beta,norm);
      auto beta = rnorm*e1
      k = 0
      while (k < maxiter - 1  and error >= tol):
  #    for k in range(0,nmax_iter-1):
          Arnoldi(Af,H,Q,k,main)
          apply_givens_rotation(H,cs,sn,k)
          #update the residual vector
          beta[k+1] = -sn[k]*beta[k]
          beta[k] = cs[k]*beta[k]
          error = abs(beta[k+1])/bnorm
          ## For testing
          #y = np.linalg.solve(H[0:k,0:k],beta[0:k]) 
          #x = x0 + np.dot(Q[:,0:k],y)
          #rt = b - Af(x)
          #rtnorm = np.linalg.norm(rt)#globalNorm(rt,main)
          if (main.mpi_rank == 0 and printnorm == 1):
            print('Outer iteration = ' + str(k_outer) + ' Iteration = ' + str(k) + '  GMRES error = ' + str(error))
            #print('Outer iteration = ' + str(k_outer) + ' Iteration = ' + str(k) + '  GMRES error = ' + str(error), ' Real norm = ' + str(rtnorm))

          k += 1
      y = np.linalg.solve(H[0:k,0:k],beta[0:k]) 
      x = x0 + np.dot(Q[:,0:k],y)
      x0[:] = x[:]
      k_outer += 1
    return x[:]
}
*/


template <typename fom_t, typename rom_state_t, typename basis_t, typename fom_state_t, typename fom_velocity_t>
class galerkinEntropySteadyRomSystem 
{
public:
  using scalar_type   = double; 
  using state_type    = rom_state_t; 
  using residual_type = rom_state_t; 
  using jacobian_type = Eigen::MatrixXd; 
  rom_state_t & romState_;
  const fom_t & appObj_;
  const basis_t & Phi_;
  mutable basis_t JPhi_;
  fom_state_t & fomState_;
  fom_velocity_t & fomVelocity_;

public:

  galerkinEntropySteadyRomSystem(const fom_t & appObj, rom_state_t & romState,const basis_t & Phi, fom_state_t & fomState, fom_velocity_t & fomVelocity) : appObj_(appObj), romState_(romState), Phi_(Phi), JPhi_(Phi) , fomState_(fomState) , fomVelocity_(fomVelocity) {};
    

  residual_type createResidual() const{
    auto result = ::pressio::ops::clone(romState_);
    ::pressio::ops::set_zero(result);
    return result;
  }

  jacobian_type createJacobian() const{
    auto K = romState_.size();
    jacobian_type reducedJacobian(K,K);
    return reducedJacobian;
  }

  void residual(const state_type & romState, residual_type & romResidual) const
  {
    ::pressio::ops::product(::pressio::nontranspose(), 1., Phi_,romState,0.,fomState_);
    appObj_.velocity(fomState_,0.,fomVelocity_);
    ::pressio::ops::product(::pressio::transpose(),1., Phi_,fomVelocity_, 0.,romResidual);    
  }
  
  void jacobian(const state_type & romState, jacobian_type & romJacobian) const
  {
    auto K = romState.size();
    /*
    Kokkos::View<scalar_type**,Kokkos::LayoutLeft> romJacobian_Kokkos("reducedMass",K,K);
    ::pressio::ops::product(::pressio::nontranspose(), 1., Phi_,romState,0.,fomState_);
    appObj_.applyJacobian(fomState_, Phi_ ,0., JPhi_);
    ::pressio::ops::product(::pressio::transpose(),::pressio::nontranspose(),1., Phi_,JPhi_, 0.,romJacobian_Kokkos);    
    for (int n=0; n < K; n++){
      for (int m=0; m< K ; m++){
        romJacobian(m,n) = romJacobian_Kokkos(m,n);
      }
    }
    */
    double eps = 1.e-4;
    auto base_residual = this->createResidual();
    auto pert_residual = this->createResidual();
    auto romStatePerturb = ::pressio::ops::clone(romState);
    this->residual(romState,base_residual);
    ::pressio::ops::set_zero(romJacobian);
    for (int i = 0; i < K; i++){
      romStatePerturb(i) += eps;
      this->residual(romStatePerturb,pert_residual);
      romStatePerturb(i) -= eps;
      for (int j = 0; j < K; j++){
        romJacobian(j,i) = ( pert_residual(j) - base_residual(j)) / eps;
      }
    } 
  }
};

template <typename fom_t, typename rom_state_t,  typename basis_t> 
class galerkinEntropySteadyRom
{
  rom_state_t & romState_;
  const fom_t & appObj_;
  const basis_t & Phi_;
public:
  galerkinEntropySteadyRom(const fom_t & appObj, rom_state_t & romState,const basis_t & Phi) : appObj_(appObj), romState_(romState), Phi_(Phi)  {};

  template <typename obs_t>
  void solveSystem(rom_state_t romState, obs_t & Obs){
    auto reducedVelocity = ::pressio::ops::clone(romState);
    auto fomState = appObj_.createVelocity();
    auto fomVelocity = appObj_.createVelocity();
    auto myGalerkinEntropySteadyRomSystem = galerkinEntropySteadyRomSystem<fom_t, rom_state_t,basis_t, decltype(fomState),decltype(fomVelocity)>(appObj_,romState,Phi_,fomState,fomVelocity);

    namespace pls  = pressio::linearsolvers;
    using tag      = pls::direct::HouseholderQR;
    using matrix_type = Eigen::MatrixXd;
    using solver_t = pls::Solver<tag, matrix_type>;
    solver_t linearSolver;

    namespace pnls = pressio::nonlinearsolvers;
    auto NonLinSolver = pnls::create_newton_raphson(myGalerkinEntropySteadyRomSystem, romState, linearSolver);
    NonLinSolver.setTolerance(1e-8);
    NonLinSolver.solve(myGalerkinEntropySteadyRomSystem,romState);
    Obs(1,0.,romState);
  }
};

}}} // end name space 

#endif
