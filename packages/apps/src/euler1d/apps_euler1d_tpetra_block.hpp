#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

namespace pressio{ namespace apps{ namespace euler1d{ namespace tpetra{


int index_map(int i1,int i2){
    return i1 + 3*i2;
}


template<typename flux_t,typename state_t>
void roeflux_kernel( flux_t  & F, const state_t & UL, const state_t & UR){
// PURPOSE: This function calculates the flux for the Euler equations
// using the Roe flux function
// INPUTS:
//    UL: conservative state vector in left cell
//   UR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)


  double gamma = 1.4;
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
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*(l3*du[0] + C1   );
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*(l3*du[1] + C1*ui + C2);
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*(l3*du[2] + C1*Hi + C2*ucp  );
  //std::cout << U(0,i) << " " << UR[0] << std::endl;

}



template<typename scalar_t>
class PressioInterface
{
  protected:
    using map_t		= Tpetra::Map<>;
    using nativeVec	= Tpetra::BlockVector<>;
    //using nativeBlockVec	= Tpetra::BlockVector<>;
 
    using go_t		= typename map_t::global_ordinal_type;
    using lo_t		= typename map_t::local_ordinal_type;
  
    using tcomm_t		= Teuchos::MpiComm<int>;
    using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
    using rcpmap_t	= Teuchos::RCP<const map_t>;

    template<typename T> using stdrcp = std::shared_ptr<T>;
    using crs_graph_type = Tpetra::CrsGraph<>;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
    rcpcomm_t comm_{};
    rcpmap_t gridMap_{};
 
    int myRank_{};
    int totRanks_{};
    lo_t NumMyElem_{};
    std::vector<go_t> myGel_{};
  
    mutable stdrcp<nativeVec> U_{}; 
    mutable stdrcp<nativeVec> U_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp2_{}; 
    stdrcp<Tpetra::Vector<>> xGrid_{}; // mesh points coordinates


//    mutable stdrcp<nativeVec> U0_{}; // initial state vector
//    stdrcp<jacobian_type> Jac_{};


  private:
    scalar_t dx_;
    //mutable state_t U_tmp_;
    //mutable state_t V_tmp_;
    //mutable state_t V_tmp2_;
    int N_cell_;
    int romSize_;

  public:
    using scalar_type	= scalar_t;
    using state_type	= nativeVec;
    using state_t = state_type;
    using velocity_type	= state_type;
    using jacobian_type	= Tpetra::CrsMatrix<>;
    using dense_matrix_type = Tpetra::BlockMultiVector<>;

    
  public:
    PressioInterface(int N_cell, scalar_t dx, int romSize, rcpcomm_t comm) : 
      N_cell_(N_cell),dx_(dx),romSize_(romSize), comm_(comm){this->setup();}
 
    //========== 
    void velocity(const state_t & U, const scalar_t & t, velocity_type & V) const {
      double U_left[3];
      double *U_left_pointer;
      double U_right[3];
      double *U_right_pointer;
      double *UL;
      double *UR;
      double *ULp;
      double *URp;
      double *V_view;

      double FL[3];
      double FR[3];


      int tag_ = 1;
      U.getLocalRowView(0,UL);
      if( myRank_ < totRanks_ - 1 ){
        U.getLocalRowView(NumMyElem_ - 1,UR);
        MPI_Send(UR, 3, MPI_DOUBLE,
  		myRank_+1, tag_, *comm_->getRawMpiComm() );
      }
      if( myRank_ > 0 ){
        MPI_Status status;
        MPI_Recv(U_left, 3, MPI_DOUBLE,
  	       myRank_-1, tag_,
  	       *comm_->getRawMpiComm(), &status);    
      }
      else{
        U_left[0] = UL[0];
        U_left[1] = -UL[1];
        U_left[2] = UL[2];
      }
      //std::cout << "test " << U_left[0]  << " " << myRank_ << std::endl;

      U_left_pointer = U_left;
      roeflux_kernel(FL,U_left_pointer,UL);
      
      for (int i=0; i < NumMyElem_ - 1 ; i++){
        U.getLocalRowView(i,UL);
        U.getLocalRowView(i+1,UR);
        roeflux_kernel(FR,UL,UR);
        V.getLocalRowView(i,V_view); 
        for (int j=0;j<3;j++){
          V_view[j] = -1./dx_*(FR[j] - FL[j]);
          FL[j] = FR[j];
        }
      }
      
      tag_ = 2;
      U.getLocalRowView(NumMyElem_ - 1,UR);
      if( myRank_ > 0){
        U.getLocalRowView(0,UL);
        MPI_Send(UL, 3, MPI_DOUBLE,
  		myRank_-1, tag_, *comm_->getRawMpiComm() );
      }
      if( myRank_ < totRanks_ - 1){
        MPI_Status status;
        MPI_Recv(U_right, 3, MPI_DOUBLE,
  	       myRank_+1, tag_,
  	       *comm_->getRawMpiComm(), &status);    
      }
      else{
        U_right[0] = UR[0];
        U_right[1] = -UR[1];
        U_right[2] = UR[2];
      }
      U_right_pointer = U_right; 
      roeflux_kernel(FR,UR,U_right_pointer);
      V.getLocalRowView(NumMyElem_ - 1,V_view); 
      for (int j=0;j<3;j++){
        V_view[j] = -1./dx_*(FR[j] - FL[j]);
      }
      
    }
    /*

    */
    
    void applyJacobian(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      double eps = 1.e-5;
      //Kokkos::deep_copy( U_tmp_, U );     
      state_t Up(*gridMap_,3);
      velocity_type V0(*gridMap_,3);
      velocity_type V_perturb(*gridMap_,3);

      velocity(U,t,V0);
      
      for (int k=0; k < romSize_; k++){
//        for (int i=0; i < N_cell_; i++){
        int i = 0;
        for (auto const & it : myGel_){
          double *Uview;
          double *Upview;
          double *Aview;
          U.getLocalRowView(i,Uview);
          Up.getLocalRowView(i,Upview);
          A.getLocalRowView(i,k,Aview);
          for (int j=0;j<3;j++){
            Upview[j] = Uview[j] + eps*Aview[j];
          }
          i++;
        }
        velocity(Up,t,V_perturb);
        i = 0;
        for (auto const & it : myGel_){
            double *V0_view;
            double *V_perturb_view;
            double *JA_view;
            V0.getLocalRowView(i,V0_view);
            V_perturb.getLocalRowView(i,V_perturb_view);
            JA.getLocalRowView(i,k,JA_view);
            for (int j=0; j < 3; j++){
              JA_view[j] = 1./eps*(V_perturb_view[j] - V0_view[j] ) ;
            }
          i++;
          }
        }
      }

    

    //========
    velocity_type createVelocity() const {
      velocity_type V(*gridMap_,3);
      return V;
    }

    dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
    {
      dense_matrix_type JA(*gridMap_, 3, A.getNumVectors() );
      return JA;
    }


    rcpmap_t getGridMap(){
      return gridMap_;
    };


   int getNumMyElem() const{
     return NumMyElem_;
   }
protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();

    // distribute cells
    gridMap_ = Teuchos::rcp(new map_t(N_cell_, 0, comm_));

    NumMyElem_ = gridMap_->getNodeNumElements();
    auto minGId = gridMap_->getMinGlobalIndex();
    myGel_.resize(NumMyElem_);
    std::iota(myGel_.begin(), myGel_.end(), minGId);


    //stateMap_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);
     // grid
    xGrid_ = std::make_shared<Tpetra::Vector<>>(gridMap_);
    auto xGridv = xGrid_->getDataNonConst();
    lo_t i=0;
    for (auto const & it : myGel_){
      xGridv[i] = dx_*it + dx_*0.5;
      i++;
    };
    //xGrid_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    // init condition
    U_ = std::make_shared<nativeVec>(*gridMap_,3);
    U_->putScalar(0.0);

    U_tmp_ = std::make_shared<nativeVec>(*gridMap_,3);
    U_tmp_->putScalar(0.0);

    V_tmp_ = std::make_shared<nativeVec>(*gridMap_,3);
    V_tmp_->putScalar(0.0);

    V_tmp2_= std::make_shared<nativeVec>(*gridMap_,3);
    V_tmp2_->putScalar(0.0);

  };

public:

  nativeVec getShockTubeIC(){
    auto xGridv = xGrid_->getData();
    double gamma = 1.4;
    double rhoL = 1;
    double pL = 1.; 
    double rhoR = 0.125;
    double pR = 0.1;
    int i = 0;
    //auto Udata = U_->getDataNonConst();
    for (auto const & it : myGel_){
      double * uVal;

      U_->getLocalRowView(i,uVal);
      if (xGridv[i] < 0.5){
        double valin[3];
        uVal[0] = rhoL;
        uVal[1]  = 0.;
        uVal[2] = pL/(gamma - 1.) ; // + 0.5*rhoU^2/rho
 
      }
      else{
        uVal[0] = rhoR;
        uVal[1]  = 0.;
        uVal[2] = pR/(gamma - 1.) ; // + 0.5*rhoU^2/rho

      }
      i++;
    }
    return *U_;
  }

};
}}}}

