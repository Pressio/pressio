#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>
//#include <Kokkos_Core.hpp>

namespace pressio{ namespace apps{ namespace euler1d{


int index_map(int i1,int i2){
    return i1 + 3*i2;
}


template<typename flux_t,typename state_t>
void roeflux_kernel(int i, flux_t  & F, const state_t & U, const int N_cell){
// PURPOSE: This function calculates the flux for the Euler equations
// using the Roe flux function
// INPUTS:
//    UL: conservative state vector in left cell
//   UR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)
  
  //Kokkos::View<double*> UR( "uR", 3);
  double UR[3];
  double UL[3];
  double gamma = 1.4;
  if (i < N_cell ){
    for (int j=0; j < 3; j++ ){
      UR[j] = U(index_map(j,i));
    }
  }
  else{
    UR[0] = U(index_map(0,N_cell-1));
    UR[1] = -U(index_map(1,N_cell-1));
    UR[2] = U(index_map(2,N_cell-1));
  }
  if (i > 0){
    for (int j=0; j<3 ;j++){
      UL[j] = U(index_map(j,i-1));
    }
  }
  else{
    UL[0] = U(index_map(0,0));
    UL[1] = -U(index_map(1,0));
    UL[2] = U(index_map(2,0));
  }
  double gmi = gamma-1.0;
  //process left state
  double rL = UL[0] + 1.e-10;
  double rhoiL = 1./rL;
  double uL = UL[1]*rhoiL;
  double unL = uL;

  double qL = pow(UL[1],2.);
  qL = pow(qL,0.5);
  qL *= rhoiL;

  double pL = (gamma-1)*(UL[2] - 0.5*rL*qL*qL);
  double rHL = UL[2] + pL;
  double HL = rHL*rhoiL;
  double cL = pow( abs(gamma*pL*rhoiL) ,0.5);
  // left flux
  double FL[ 3 ];
  FL[0] = rL*unL;
  FL[1] = UL[1]*unL + pL;
  FL[2] = rHL*unL;

  // process right state
  double rR = UR[0];
  double rhoiR = 1./rR;
  double uR = UR[1]*rhoiR;

  double unR =  uR;

  double qR = pow(UR[1],2);
  qR = pow(qR,0.5);
  qR *= rhoiR;

  double pR = (gamma-1)*(UR[2]- 0.5*rR*qR*qR);
  double rHR = UR[2] + pR;
  double HR = rHR*rhoiR;
  double cR = pow(abs(gamma*pR*rhoiR),0.5);
  // right flux
  double FR[3];
  FR[0] = rR*unR;
  FR[1] = UR[1]*unR + pR;
  FR[2] = rHR*unR;
  // difference in states
  double du[3]; 
  du[0] = UR[0] - UL[0];
  du[1] = UR[1] - UL[1];
  du[2] = UR[2] - UL[2];

  // Roe average
  double di     = pow((rR*rhoiL),0.5);
  double d1     = 1.0/(1.0+di);

  double ui     = (di*uR + uL)*d1;
  double Hi     = (di*HR + HL)*d1;

  double af     = 0.5*(ui*ui);
  double ucp    = ui;
  double c2     = gmi*(Hi - af);
  double ci     = pow(abs(c2),0.5);
  double ci1    = 1.0/ci;


  // eigenvalues
  double l[ 3 ];
  l[0] = ucp+ci;
  l[1] = ucp-ci;
  l[2] = ucp;

  // entropy fix
  double epsilon = ci*.1;
  double labs[ 3]; 
  for (int k=0 ;  k<3 ;  k++){ labs[k] =  abs(l[k]);}
  for (int k = 0; k < 3; k++){ 
    if (labs[k] < epsilon){
    l[k] =  (epsilon + l[k]*l[k])/(2.*epsilon);
    }
  }

  for (int k=0 ;  k<3 ;  k++){ l[k] =  abs(l[k]);}

  double l3 = l[2];
  //average and half-difference of 1st and 2nd eigs
  double s1    = 0.5*(l[0] + l[1]);
  double s2    = 0.5*(l[0] - l[1]);

  //left eigenvector product generators (see Theory guide)
  double G1    = gmi*(af*du[0] - ui*du[1] + du[2]);
  double G2    = -ucp*du[0]+du[1];
  //required functions of G1 and G2 (again, see Theory guide)
  double C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
  double C2    = G1*s2*ci1          + G2*(s1-l3);
  
  //# flux assembly
  F(index_map(0,i))    = 0.5*(FL[0]+FR[0])-0.5*(l3*du[0] + C1   );
  F(index_map(1,i))    = 0.5*(FL[1]+FR[1])-0.5*(l3*du[1] + C1*ui + C2);
  F(index_map(2,i))    = 0.5*(FL[2]+FR[2])-0.5*(l3*du[2] + C1*Hi + C2*ucp  );
  //std::cout << U(0,i) << " " << UR[0] << std::endl;
}

template <typename velocity_t, typename flux_t, typename state_t>
void computeVelocity(velocity_t & velocity, flux_t  & F, const state_t & U, int N_cell, double dx){

  Kokkos::parallel_for( "roe_flux", N_cell + 1, KOKKOS_LAMBDA ( int i ) {
    roeflux_kernel(i,F,U,N_cell);
  });

  Kokkos::parallel_for( "velocity", N_cell , KOKKOS_LAMBDA ( int i ) {
      for (int j=0;j<3;j++){
        velocity(index_map(j,i)) = -1./dx*(F(index_map(j,i+1)) - F(index_map(j,i)));
      }
  });
}


template<typename state_t, typename velocity_t, typename flux_t, typename scalar_t>
class PressioInterface
{
  private:
    scalar_t dx_;
    mutable state_t U_tmp_;
    mutable state_t V_tmp_;
    mutable state_t V_tmp2_;
    int N_cell_;
    int N_hyper_ = 153;
    int romSize_;
    typedef Kokkos::LayoutLeft   Layout;
    using matrix_t = Kokkos::View<double**,Layout>;
    
  public:
    flux_t flux_;
    using scalar_type = double;
    using state_type = state_t;
    using velocity_type = state_t;
    using jacobian_type = matrix_t;
    using dense_matrix_type = matrix_t;
    int index_mappings [153] = { 0,    1,    2,   30,   31,   32,   60,   61,   62,   90,   91,
         92,  120,  121,  122,  150,  151,  152,  180,  181,  182,  210,
        211,  212,  240,  241,  242,  270,  271,  272,  300,  301,  302,
        330,  331,  332,  360,  361,  362,  390,  391,  392,  420,  421,
        422,  450,  451,  452,  480,  481,  482,  510,  511,  512,  540,
        541,  542,  570,  571,  572,  600,  601,  602,  630,  631,  632,
        660,  661,  662,  690,  691,  692,  720,  721,  722,  750,  751,
        752,  780,  781,  782,  810,  811,  812,  840,  841,  842,  870,
        871,  872,  900,  901,  902,  930,  931,  932,  960,  961,  962,
        990,  991,  992, 1020, 1021, 1022, 1050, 1051, 1052, 1080, 1081,
       1082, 1110, 1111, 1112, 1140, 1141, 1142, 1170, 1171, 1172, 1200,
       1201, 1202, 1230, 1231, 1232, 1260, 1261, 1262, 1290, 1291, 1292,
       1320, 1321, 1322, 1350, 1351, 1352, 1380, 1381, 1382, 1410, 1411,
       1412, 1440, 1441, 1442, 1470, 1471, 1472, 1497, 1498, 1499};


    

    PressioInterface(int N_cell, scalar_t dx, int romSize) : 
      U_tmp_("Utmp",N_cell*3),V_tmp_("Vtmp",N_cell*3),V_tmp2_("Vtmp2",N_cell*3),N_cell_(N_cell) , dx_(dx), flux_("Flux",(N_cell_+1)*3), romSize_(romSize){}
 

    //========== 
    void velocity(const state_t & U, const scalar_t & t, velocity_t & V) const {
      state_t flux("f",3*(N_cell_+1));
      computeVelocity(V,flux_,U,N_cell_,dx_);
    }

    velocity_t createVelocity() const {
      state_t V("VelocityF",3*N_cell_);
      return V;
    }


    void applyJacobian(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      //Kokkos::deep_copy( U_tmp_, U );     
      velocity(U,t,V_tmp2_);
      for (int k=0; k < romSize_; k++){
        for (int i=0; i < 3; i++){
          for (int j=0; j < N_cell_; j++){
            U_tmp_(index_map(i,j)) = U(index_map(i,j)) + eps*A(index_map(i,j),k);}}

        velocity(U_tmp_,t,V_tmp_);

        for (int i=0; i < 3; i++){
          for (int j=0; j < N_cell_; j++){
            JA(index_map(i,j),k) = 1./eps*(V_tmp_(index_map(i,j)) - V_tmp2_(index_map(i,j)) ) ;}}
      }
    }

    dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A)const{
      dense_matrix_type JA("JA",3*N_cell_,romSize_);
      return JA;
    }


    //========

    //===== hyper reduction ///
    int hyper_map(const int i) const {
      return index_mappings[i];
    }
 
    state_t get_reduced_state(const state_t & u){
      state_t u_hyper("u_h",N_hyper_);
      for (int i = 0; i < N_hyper_; i++){
        u_hyper(i) = u( hyper_map(i));
      }
      return u_hyper;
    } 
    
    void velocity2(const state_t & U, const scalar_t & t, velocity_t & V) const {
      state_t flux("f",N_hyper_);
      computeVelocity(V_tmp_,flux_,U,N_cell_,dx_);
      for (int i = 0; i < N_hyper_; i++){
        V(i) = V_tmp_( hyper_map(i));
      }
    }

    velocity_t createVelocity2() const {
      state_t V("VelocityF",N_hyper_);
      return V;
    }


    void applyJacobian2(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      state_t v1_hyper("v1h",N_hyper_);
      state_t v2_hyper("v2h",N_hyper_);
      Kokkos::deep_copy(U_tmp_, U);
      velocity(U,t,v1_hyper);
      for (int k=0; k < romSize_; k++){
        for (int i=0; i < N_hyper_; i++){
            U_tmp_( hyper_map(i) ) = U( hyper_map(i) ) + eps*A(i,k);}

        velocity(U_tmp_,t,v2_hyper);

        for (int i=0; i < N_hyper_; i++){
          JA(i,k) = 1./eps*(v2_hyper(i) - v1_hyper(i) ) ;
          }
      }
    }

    dense_matrix_type createApplyJacobianResult2(const dense_matrix_type & A)const{
      dense_matrix_type JA("JA",N_hyper_,romSize_);
      return JA;
    }


    //===========



};
}}}

