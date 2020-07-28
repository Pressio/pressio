#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

namespace pressio{ namespace apps{ namespace euler1d{ namespace tpetra{ namespace hyper{




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

  private:
    scalar_t dx_;
    //mutable state_t U_tmp_;
    //mutable state_t V_tmp_;
    //mutable state_t V_tmp2_;
    int N_cell_;
    int romSize_;
    const int N_sample_mesh_plus_stencil_ = 131;
    const int N_sample_mesh_ = 50;


  protected:
    Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
    rcpcomm_t comm_{};
    rcpmap_t hyperMap_{};
    rcpmap_t hyperMapWithStencil_{};

    int myRank_{};
    int totRanks_{};
    lo_t NumMyElem_{};
    lo_t NumMyElemHyper_{};

    std::vector<go_t> myGel_{};
    std::vector<go_t> myGelHyper_{};
 

    mutable stdrcp<nativeVec> Us_{}; 


    mutable stdrcp<nativeVec> U_{}; 
    mutable stdrcp<nativeVec> U_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp_{}; 
    mutable stdrcp<nativeVec> V_tmp2_{}; 
    stdrcp<Tpetra::Vector<>> xGrid_{}; // mesh points coordinates
    Kokkos::View<go_t*> sampleIndices_; 



  public:
    using scalar_type	= scalar_t;
    using state_type	= nativeVec;
    using state_t = state_type;
    using velocity_type	= state_type;
    using jacobian_type	= Tpetra::CrsMatrix<>;
    using dense_matrix_type = Tpetra::BlockMultiVector<>;

    
  public:
    PressioInterface(int N_cell, scalar_t dx, int romSize, rcpcomm_t comm) : 
      sampleIndices_("siv",N_sample_mesh_), N_cell_(N_cell),dx_(dx),romSize_(romSize), comm_(comm){this->setup();}
 
    //========== 
    void velocity(const state_t & U, const scalar_t & t, velocity_type & V) const {
      // U is on the sample mesh
      //  note that U is ordered on the sample mesh
      // V is not
      double U_left[3];
      double *U_left_pointer = nullptr;
      double U_right[3];
      double *U_right_pointer = nullptr;
      double *UL = nullptr;
      double *UR = nullptr;
      double *V_view = nullptr;

      double FL[3];
      double FR[3];

      // First grid point is on the sample mesh by design
      U.getLocalRowView(0,UL);
      U_left[0] = UL[0];
      U_left[1] = -UL[1];
      U_left[2] = UL[2];
      U_left_pointer = U_left;
      roeflux_kernel(FL,U_left_pointer,UL);
      U.getLocalRowView(1,UR);
      roeflux_kernel(FR,UL,UR);
      V.getLocalRowView(0,V_view);      
      for (int j=0;j<3;j++){
        V_view[j] = -1./dx_*(FR[j] - FL[j]);
      }
      //================= 
      for (int i=1; i < N_sample_mesh_ - 1 ; i++){
        U.getLocalRowView(sampleIndices_(i-1),UL);
        U.getLocalRowView(sampleIndices_(i),UR);
        roeflux_kernel(FL,UL,UR);
        U.getLocalRowView(sampleIndices_(i),UL);
        U.getLocalRowView(sampleIndices_(i+1),UR);
        roeflux_kernel(FR,UL,UR);
        V.getLocalRowView(i,V_view); 
        for (int j=0;j<3;j++){
          V_view[j] = -1./dx_*(FR[j] - FL[j]);
        }
      }

      U.getLocalRowView(sampleIndices_(N_sample_mesh_-2),UL);
      U.getLocalRowView(sampleIndices_(N_sample_mesh_-1),UR);
      U_right[0] = UR[0];
      U_right[1] = -UR[1];
      U_right[2] = UR[2];
      U_right_pointer = U_right; 
      roeflux_kernel(FR,UR,U_right_pointer);
      V.getLocalRowView(N_sample_mesh_ - 1,V_view); 
      for (int j=0;j<3;j++){
        V_view[j] = -1./dx_*(FR[j] - FL[j]);
      }
    }

    
    void applyJacobian(const state_t &U, const dense_matrix_type & A,scalar_t t, dense_matrix_type &JA)const{
      // A will be on the sample mesh with stencil
      double eps = 1.e-5;
      state_t Up(*hyperMapWithStencil_,3);
      velocity_type V0(*hyperMap_,3);
      velocity_type V_perturb(*hyperMap_,3);
      Up.putScalar(0.);
      V_perturb.putScalar(0.);

      velocity(U,t,V0);
      double *Av;
      A.getLocalRowView(sampleIndices_(0),1,Av);
      for (int k=0; k < romSize_; k++){

        for (int i=0; i < N_sample_mesh_  ; i++){

        double *Uview = nullptr;
        double *Upview = nullptr;
        double *Aview = nullptr;

          U.getLocalRowView(sampleIndices_(i),Uview);
          Up.getLocalRowView(sampleIndices_(i),Upview);
          A.getLocalRowView(sampleIndices_(i),k,Aview);
          for (int j=0;j<3;j++){
            Upview[j] = Uview[j] + eps*Aview[j];
          }
        }
        velocity(Up,t,V_perturb);

        for (int i=0; i < N_sample_mesh_  ; i++){
          double *V0_view = nullptr;
          double *V_perturb_view = nullptr;
          double *JA_view = nullptr;
          V0.getLocalRowView(i,V0_view);
          V_perturb.getLocalRowView(i,V_perturb_view);
          JA.getLocalRowView(i,k,JA_view);
          for (int j=0; j < 3; j++){
            JA_view[j] = 1./eps*(V_perturb_view[j] - V0_view[j] ) ;
          }
        }

      }
    }


    //========
    velocity_type createVelocity() const {
      velocity_type V(*hyperMap_,3);
      return V;
    }

    dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
    {
      dense_matrix_type JA(*hyperMap_, 3, A.getNumVectors() );
      return JA;
    }


    rcpmap_t getHyperMapWithStencil(){
      return hyperMapWithStencil_;
    };

protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();

    // distribute cells
    //gridMap_ = Teuchos::rcp(new map_t(N_cell_, 0, comm_));



    //stateMap_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);
    //xGrid_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    // init condition


    std::ifstream sample_mesh_ind_file("sample_mesh_relative_indices.txt");
  
    for (int i =0; i < N_sample_mesh_ ; i++){
      sample_mesh_ind_file >> sampleIndices_(i);
    }
    sample_mesh_ind_file.close();
  
    //Teuchos::ArrayView<int> sampleIndicesView( sampleIndices_);
  
    hyperMap_ = Teuchos::rcp(new map_t(N_sample_mesh_,sampleIndices_,0, comm_));
    hyperMapWithStencil_ = Teuchos::rcp(new map_t(N_sample_mesh_plus_stencil_,0, comm_));
    //hyperMap_ = Tpetra::createNonContigMap<lo_t,go_t>(sampleIndices_,comm_);
    Us_ = std::make_shared<nativeVec>(*hyperMapWithStencil_,3);
    U_ = std::make_shared<nativeVec>(*hyperMap_,3);
    U_->putScalar(1.);
    Us_->putScalar(2.);

    double vals[3]; 
    vals[0] = 11;
    vals[1] = 12;
    vals[2] = 13;
    double * valsp = vals;
    U_->replaceLocalValues(2,valsp);
    NumMyElem_ = hyperMapWithStencil_->getNodeNumElements();
    auto minGId = hyperMapWithStencil_->getMinGlobalIndex();
    myGel_.resize(NumMyElem_);
    std::iota(myGel_.begin(), myGel_.end(), minGId);


    NumMyElemHyper_ = hyperMap_->getNodeNumElements();
    auto minGIdHyper = hyperMap_->getMinGlobalIndex();
    myGelHyper_.resize(NumMyElemHyper_);
    std::iota(myGelHyper_.begin(), myGelHyper_.end(), minGIdHyper);
    for (auto it : myGelHyper_){
       double * test;
       U_->getGlobalRowView(sampleIndices_(it),test);
     }




    };

public:

  //template<typename x_typ
  nativeVec getShockTubeIC(){
    std::ifstream U0_file("U0_SampleMeshPlusStencil.txt");
    for (int i=0; i < N_sample_mesh_plus_stencil_; i++){
      double * Uview;
      Us_->getLocalRowView(i,Uview);
      for (int j=0; j < 3; j++){
        U0_file >> Uview[j];
      }  
    }
    U0_file.close();
    return *Us_;
  }
};
}}}}}

