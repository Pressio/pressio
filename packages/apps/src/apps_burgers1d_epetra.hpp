
#ifndef APPS_BURGERS1D_EPETRA_HPP_
#define APPS_BURGERS1D_EPETRA_HPP_

#include "apps_ConfigDefs.hpp"
#include <memory>
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"


namespace apps{

class Burgers1dEpetra{
private:
  using nativeVec = Epetra_Vector;
  template<typename T>
  using rcp = std::shared_ptr<T>;

public:
  using scalar_type = double;
  using state_type = Epetra_Vector;
  using space_residual_type = Epetra_Vector;
  using jacobian_type = Epetra_CrsMatrix;

public:  
  Burgers1dEpetra(std::vector<scalar_type> params,
		  int Ncell, Epetra_MpiComm * comm)
    : mu_(params), Ncell_(Ncell), comm_(comm){}

  ~Burgers1dEpetra() = default; 
 
  Epetra_Map const & getDataMap(){
    return *dataMap_;
  };
   
  void setup(){
    // distribute cells 
    dataMap_ = std::make_shared<Epetra_Map>(Ncell_,0,*comm_);

    NumMyElem_ = dataMap_->NumMyElements();
    myGel_.resize(NumMyElem_);
    dataMap_->MyGlobalElements(myGel_.data());

    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;
    
    // grid 
    xGrid_ = std::make_shared<nativeVec>(*dataMap_);
    int i=0;
    for (auto const & it : myGel_){
      (*xGrid_)[i] = dx_*it + dx_*0.5;
      i++;
    };
    //xGrid_->Print(std::cout);
    // init condition
    U0_ = std::make_shared<nativeVec>(*dataMap_);
    U_ = std::make_shared<nativeVec>(*dataMap_);
    U0_->PutScalar(1.0);
    U_->PutScalar(1.0);
    myRank_ =  comm_->MyPID();
    totRanks_ =  comm_->NumProc();
    
  };

  state_type const & getInitialState(){
    return *U0_;
  };

  space_residual_type const & getInitialResidual(){
    r0_ = std::make_shared<nativeVec>(*dataMap_);
    this->residual(*U0_, *r0_, 0.0);
    return *r0_;
  };
  
  void residual(const state_type & u,
		space_residual_type & rhs,
		const scalar_type /* t */)
  {
    double valueFromLeft = 0.0;
    constexpr int tag_ = 1;
    if( myRank_ < comm_->NumProc()-1 ){
      MPI_Send( &u[NumMyElem_-1],1, MPI_DOUBLE,
	        myRank_+1, tag_, comm_->Comm() );
    }    
    if( myRank_ > 0 ){
      MPI_Status status;
      MPI_Recv(&valueFromLeft,1, MPI_DOUBLE,myRank_-1,
	       tag_, comm_->Comm(), &status);
    }
    
    int i=0;
    scalar_type uim1;
    for (auto const & it : myGel_){
      uim1 = valueFromLeft;
      if (it==0)
	uim1 = mu_[0]; // left boundary condition
      if (i>0)
	uim1 = u[i-1];

      rhs[i] = ( 0.5*(uim1*uim1 - u[i]*u[i]) )/dx_;
      i++;
    }

    for (i=0; i<NumMyElem_; ++i){
      rhs[i] += mu_[1]*exp(mu_[2] * (*xGrid_)[i]);
    }  
  }//end residual

  jacobian_type const & getInitialJacobian(){
    j0_ = std::make_shared<Epetra_CrsMatrix>(Copy, *dataMap_, nonZrPerRow_);
    this->jacobian(*U0_, *j0_, 0.0);
    j0_->FillComplete();
    return *j0_;
  };

  
  void jacobian(const state_type & u,
		jacobian_type & jac,
		const scalar_type /*t*/)
  {

    // to populate the jacobian each process needs the last grid
    // point solution from the previous process
    double buffin = 0.0;
    MPI_Comm mpiComm = comm_->Comm();
    if (myRank_ < totRanks_-1){
      double tosend = u[NumMyElem_-1];
      MPI_Send(&tosend, 1, MPI_DOUBLE, myRank_+1, 1, mpiComm);
    }
    if (myRank_ >= 1){
      MPI_Status st;
      MPI_Recv(&buffin, 1, MPI_DOUBLE, myRank_-1, 1, mpiComm, &st);
    }

    std::vector<int> Indices {0, 0};
    std::vector<double> Values {0., 0.};
    for (int i=0; i<NumMyElem_; i++)
    {
       int thisGID = myGel_[i]; // global ID
       if (thisGID==0){
       	 Indices[0] = 0;
       	 Values[0] = -dxInv_;// * u[0]*u[0];
       	 jac.InsertGlobalValues(thisGID, 1, Values.data(), Indices.data() );
       }
       else{
       	 Indices[0] = thisGID-1;
       	 Indices[1] = thisGID;
    	 if (i==0)
    	   Values[0] = dxInv_ * buffin*buffin;
    	 if (i>0)
    	   Values[0] = dxInv_ * u[i-1]*u[i-1];

    	 Values[1] = -Values[0];	 
       	 jac.InsertGlobalValues(thisGID, 2, Values.data(), Indices.data() );
       }
    }
    if (!jac.Filled())
      jac.FillComplete();
  }//end jacobian
  
private:  
  std::vector<scalar_type> mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain 
  const scalar_type xR_ = 100.0; // right side of domain
  int Ncell_; // # of cells
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size
  const int nonZrPerRow_ = 2;
  // mesh points coordinates
  rcp<nativeVec> xGrid_; 

  Epetra_MpiComm * comm_;
  rcp<Epetra_Map> dataMap_;
  int myRank_;
  int totRanks_;
  int NumMyElem_;
  std::vector<int> myGel_;  
  
  rcp<nativeVec> U_; // state vector
  rcp<nativeVec> U0_; // initial state vector
  rcp<nativeVec> r0_; // initial space residual
  rcp<Epetra_CrsMatrix> j0_; // initial jacobian
};

}//end namespace apps
#endif 
