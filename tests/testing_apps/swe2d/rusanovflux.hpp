
#ifndef RUSANOV_FLUX_HPP_
#define RUSANOV_FLUX_HPP_

#include <math.h>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

/* Rusanov flux scheme that takes the entire state vector as an input
F \in R^3, flux
q \in R^N, state
isL \in NaturalNumbers, left index
isR \in NaturalNumbers, right index
n \in R^2, normalv ector
g \in R, gravity
*/
template<typename flux_t,typename state_t, typename normal_t, typename scalar_t>
void rusanovFluxFullStateIn( flux_t  & F, const state_t & q,
		     const int isL, const int isR, const normal_t n, const scalar_t g)
{
// Computes the flux for the shallow water equations  
  const scalar_t es = 1.e-30;
  const scalar_t hL = q(isL + 0);
  const scalar_t uL = q(isL + 1)/(hL + es);
  const scalar_t vL = q(isL + 2)/(hL + es);
  const scalar_t unL = uL*n[0] + vL*n[1];

  const scalar_t pL = 0.5*g*std::pow(hL,2.);
  scalar_t FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = q(isL + 1)*unL + pL*n[0];
  FL[2] = q(isL + 2)*unL + pL*n[1];

  const scalar_t hR = q(isR + 0);
  const scalar_t uR = q(isR + 1)/(hR + es);
  const scalar_t vR = q(isR + 2)/(hR + es);
  const scalar_t unR = uR*n[0] + vR*n[1];
  const scalar_t pR = 0.5*g*std::pow(hR,2.);
  scalar_t FR[3];
  FR[0] = hR*unR;
  FR[1] = q(isR + 1)*unR + pR*n[0];
  FR[2] = q(isR + 2)*unR + pR*n[1];

  // rho average
  const scalar_t hm = 0.5*(hL + hR);
  const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));
  // flux assembly
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(q(isR + 0) - q(isL + 0));
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(q(isR + 1) - q(isL + 1));
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(q(isR + 2) - q(isL + 2));
}

/* Same as above, but takes in the left state and right state as inputs
F \in R^3, flux
qL \in R^3, left state
qR \in R^3, right state
n \in R^2, normalv ector
g \in R, gravity
*/
template<typename flux_t,typename state_t, typename normal_t, typename scalar_t>
void rusanovFluxLeftRightStateIn( flux_t  & F, const state_t & qL,
		     const state_t & qR, const normal_t n, const scalar_t g)
{
// Computes the flux for the shallow water equations  
  const scalar_t es = 1.e-30;
  const scalar_t hL = qL[0];
  const scalar_t uL = qL[1]/(hL + es);
  const scalar_t vL = qL[2]/(hL + es);
  const scalar_t unL = uL*n[0] + vL*n[1];

  const scalar_t pL = 0.5*g*std::pow(hL,2.);
  scalar_t FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = qL[1]*unL + pL*n[0];
  FL[2] = qL[2]*unL + pL*n[1];

  const scalar_t hR = qR[0];
  const scalar_t uR = qR[1]/(hR + es);
  const scalar_t vR = qR[2]/(hR + es);
  const scalar_t unR = uR*n[0] + vR*n[1];
  const scalar_t pR = 0.5*g*std::pow(hR,2.);
  scalar_t FR[3];
  FR[0] = hR*unR;
  FR[1] = qR[1]*unR + pR*n[0];
  FR[2] = qR[2]*unR + pR*n[1];

  // rho average
  const scalar_t hm = 0.5*(hL + hR);
  const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));
  // flux assembly
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(qR[0] - qL[0]);
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(qR[1] - qL[1]);
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(qR[2] - qL[2]);

}



/* Rusanov flux Jacobian that takes the entire state vector as an input
JL \in R^{3x3}, left jacobian
JR \in R^{3x3}, right jacobian
q \in R^N, state
isL \in NaturalNumbers, left index
isR \in NaturalNumbers, right index
n \in R^2, normalv ector
g \in R, gravity
*/
template<typename jac_t,typename state_t, typename normal_t, typename scalar_t>
void rusanovFluxJacobianFullStateIn( jac_t  & JL , jac_t & JR, const state_t & q,
		     const int gidL, const int gidR, const normal_t n, const scalar_t g){

// Computes the flux Jacobian for the shallow water equations  
  const scalar_t es = 1.e-30;
  const scalar_t hL = q(gidL + 0);
  const scalar_t uL = q(gidL + 1)/(hL + es);
  const scalar_t vL = q(gidL + 2)/(hL + es);
  const scalar_t unL = uL*n[0] + vL*n[1];

  const scalar_t hR = q(gidR + 0);
  const scalar_t uR = q(gidR + 1)/(hR + es);
  const scalar_t vR = q(gidR + 2)/(hR + es);
  const scalar_t unR = uR*n[0] + vR*n[1];

  // rho average
  const scalar_t hm = 0.5*(hL + hR);
  const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));


  const scalar_t termL = (n[0]*q(gidL + 1) + n[1]*q(gidL + 2)) / std::pow(q(gidL + 0),2.);

  const scalar_t termR = (n[0]*q(gidR + 1) + n[1]*q(gidR + 2)) / std::pow(q(gidR + 0),2.);

  const scalar_t hL_sqrt = std::pow(hL,0.5);
  const scalar_t hR_sqrt = std::pow(hR,0.5);

  const scalar_t hsqrt_un = hL_sqrt*unL + hR_sqrt*unR + es;
  scalar_t dsmaxL[3];
  scalar_t dsmaxR[3];
  dsmaxL[0] = - std::abs( hsqrt_un ) / (2.*hL_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unL / hL_sqrt - hL_sqrt*termL )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un  ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxL[1] = n[0]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxL[2] = n[1]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );


  dsmaxR[0] = - std::abs( hsqrt_un ) / (2.*hR_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unR / hR_sqrt - hR_sqrt*termR )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxR[1] = n[0]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxR[2] = n[1]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  // jacobian w.r.p to the left state
  JL[0][0] = -0.5*dsmaxL[0]*(q(gidR + 0) - q(gidL + 0)) + 0.5*(n[0]*uL + n[1]*vL - q(gidL + 0)*termL)  + 0.5*smax;
  JL[0][1] = 0.5*n[0] - 0.5*dsmaxL[1]*(q(gidR + 0) - q(gidL + 0));
  JL[0][2] = 0.5*n[1] - 0.5*dsmaxL[2]*(q(gidR + 0) - q(gidL + 0));

  JL[1][0] = 0.5*(g*n[0]*q(gidL + 0) - q(gidL + 1)*termL) - 0.5*dsmaxL[0]*(q(gidR + 1) - q(gidL + 1));
  JL[1][1] = n[0]*uL  + 0.5*n[1]*vL + 0.5*smax -0.5*dsmaxL[1]*(q(gidR + 1) - q(gidL + 1));
  JL[1][2] = 0.5*n[1]*uL  -0.5*dsmaxL[2]*(q(gidR + 1) - q(gidL + 1));

  JL[2][0] = 0.5*(g*n[1]*q(gidL + 0) - q(gidL + 2)*termL) - 0.5*dsmaxL[0]*(q(gidR + 2) - q(gidL + 2));
  JL[2][1] = 0.5*n[0]*vL  -0.5*dsmaxL[1]*(q(gidR + 2) - q(gidL + 2));
  JL[2][2] = n[1]*vL + 0.5*n[0]*uL + 0.5*smax -0.5*dsmaxL[2]*(q(gidR + 2) - q(gidL + 2));

   // jacobian w.r.p to the right state
  JR[0][0] = -0.5*dsmaxR[0]*(q(gidR + 0) - q(gidL + 0)) + 0.5*(n[0]*uR + n[1]*vR - q(gidR + 0)*termR)  - 0.5*smax;
  JR[0][1] = 0.5*n[0] - 0.5*dsmaxR[1]*(q(gidR + 0) - q(gidL + 0));
  JR[0][2] = 0.5*n[1] - 0.5*dsmaxR[2]*(q(gidR + 0) - q(gidL + 0));

  JR[1][0] = 0.5*(g*n[0]*q(gidR + 0) - q(gidR + 1)*termR) - 0.5*dsmaxR[0]*(q(gidR + 1) - q(gidL + 1));
  JR[1][1] = n[0]*uR + 0.5*n[1]*vR - 0.5*smax - 0.5*dsmaxR[1]*(q(gidR + 1) - q(gidL + 1));
  JR[1][2] = 0.5*n[1]*uR -0.5*dsmaxR[2]*(q(gidR + 1) - q(gidL + 1));

  JR[2][0] = 0.5*(g*n[1]*q(gidR + 0)  - q(gidR + 2)*termR) - 0.5*dsmaxR[0]*(q(gidR + 2) - q(gidL + 2));
  JR[2][1] = 0.5*n[0]*vR - 0.5*dsmaxR[1]*(q(gidR + 2) - q(gidL + 2));
  JR[2][2] = n[1]*vR + 0.5*n[0]*uR - 0.5*smax - 0.5*dsmaxR[2]*(q(gidR + 2) - q(gidL + 2));

}



/* Rusanov flux Jacobian that takes the left and right state as an input
JL \in R^{3x3}, left jacobian
JR \in R^{3x3}, right jacobian
qL \in R^3, left state
qR \in R^3, right state
n \in R^2, normalv ector
g \in R, gravity
*/
template<typename jac_t,typename state_t, typename normal_t, typename scalar_t>
void rusanovFluxJacobianLeftRightStateIn( jac_t  & JL , jac_t & JR, const state_t & qL,
		     const state_t & qR, const normal_t n,const scalar_t g){

// Computes the flux Jacobian for the shallow water equations  
  const scalar_t es = 1.e-30;

  const scalar_t hL = qL[0];
  const scalar_t uL = qL[1]/(hL + es);
  const scalar_t vL = qL[2]/(hL + es);
  const scalar_t unL = uL*n[0] + vL*n[1];

  const scalar_t hR = qR[0];
  const scalar_t uR = qR[1]/(hR + es);
  const scalar_t vR = qR[2]/(hR + es);
  const scalar_t unR = uR*n[0] + vR*n[1];

  // rho average
  const scalar_t hm = 0.5*(hL + hR);
  const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));


  const scalar_t termL = (n[0]*qL[1] + n[1]*qL[2]) / std::pow(qL[0],2.);

  const scalar_t termR = (n[0]*qR[1] + n[1]*qR[2]) / std::pow(qR[0],2.);

  const scalar_t hL_sqrt = std::pow(hL,0.5);
  const scalar_t hR_sqrt = std::pow(hR,0.5);

  const scalar_t hsqrt_un = hL_sqrt*unL + hR_sqrt*unR + es;
  scalar_t dsmaxL[3];
  scalar_t dsmaxR[3];
  dsmaxL[0] = - std::abs( hsqrt_un ) / (2.*hL_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unL / hL_sqrt - hL_sqrt*termL )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un  ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxL[1] = n[0]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxL[2] = n[1]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );


  dsmaxR[0] = - std::abs( hsqrt_un ) / (2.*hR_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unR / hR_sqrt - hR_sqrt*termR )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxR[1] = n[0]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxR[2] = n[1]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  // jacobian w.r.p to the left state
  JL[0][0] = -0.5*dsmaxL[0]*(qR[0] - qL[0]) + 0.5*(n[0]*uL + n[1]*vL - qL[0]*termL)  + 0.5*smax;
  JL[0][1] = 0.5*n[0] - 0.5*dsmaxL[1]*(qR[0] - qL[0]);
  JL[0][2] = 0.5*n[1] - 0.5*dsmaxL[2]*(qR[0] - qL[0]);

  JL[1][0] = 0.5*(g*n[0]*qL[0] - qL[1]*termL) - 0.5*dsmaxL[0]*(qR[1] - qL[1]);
  JL[1][1] = n[0]*uL  + 0.5*n[1]*vL + 0.5*smax -0.5*dsmaxL[1]*(qR[1] - qL[1]);
  JL[1][2] = 0.5*n[1]*uL  -0.5*dsmaxL[2]*(qR[1] - qL[1]);

  JL[2][0] = 0.5*(g*n[1]*qL[0] - qL[2]*termL) - 0.5*dsmaxL[0]*(qR[2] - qL[2]);
  JL[2][1] = 0.5*n[0]*vL  -0.5*dsmaxL[1]*(qR[2] - qL[2]);
  JL[2][2] = n[1]*vL + 0.5*n[0]*uL + 0.5*smax -0.5*dsmaxL[2]*(qR[2] - qL[2]);

   // jacobian w.r.p to the right state
  JR[0][0] = -0.5*dsmaxR[0]*(qR[0] - qL[0]) + 0.5*(n[0]*uR + n[1]*vR - qR[0]*termR)  - 0.5*smax;
  JR[0][1] = 0.5*n[0] - 0.5*dsmaxR[1]*(qR[0] - qL[0]);
  JR[0][2] = 0.5*n[1] - 0.5*dsmaxR[2]*(qR[0] - qL[0]);

  JR[1][0] = 0.5*(g*n[0]*qR[0] - qR[1]*termR) - 0.5*dsmaxR[0]*(qR[1] - qL[1]);
  JR[1][1] = n[0]*uR + 0.5*n[1]*vR - 0.5*smax - 0.5*dsmaxR[1]*(qR[1] - qL[1]);
  JR[1][2] = 0.5*n[1]*uR -0.5*dsmaxR[2]*(qR[1] - qL[1]);

  JR[2][0] = 0.5*(g*n[1]*qR[0]  - qR[2]*termR) - 0.5*dsmaxR[0]*(qR[2] - qL[2]);
  JR[2][1] = 0.5*n[0]*vR - 0.5*dsmaxR[1]*(qR[2] - qL[2]);
  JR[2][2] = n[1]*vR + 0.5*n[0]*uR - 0.5*smax - 0.5*dsmaxR[2]*(qR[2] - qL[2]);
}

#endif
