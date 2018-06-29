
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

class burgers1dEpetra
{
private:
  using nativeVec = Epetra_Vector;

  template<typename T>
  using rcp = std::shared_ptr<T>;

public:
  using scalar_type = double;
  using state_type = nativeVec;
  using jacobian_type = Epetra_CrsMatrix;

public:  
  burgers1dEpetra(std::vector<scalar_type> params,
		  Epetra_MpiComm * comm,
		  int Ncell=1000)
    : mu_(params), comm_(comm), Ncell_(Ncell)
  {}

  ~burgers1dEpetra() = default; 

  void setup()
  {
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
    xGrid_->Print(std::cout);
    // init condition
    U0_ = std::make_shared<nativeVec>(*dataMap_);
    U_ = std::make_shared<nativeVec>(*dataMap_);
    U0_->PutScalar(1.0);
    U_->PutScalar(1.0);

    myRank_ =  comm_->MyPID();
  };

  state_type getInitialState(){
    return *U0_;
  };
  
  void residual(const state_type & u,
		state_type & rhs,
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
    scalar_type uim1, ui;
    for (auto const & it : myGel_)
    {
      uim1 = valueFromLeft;
      if (it==0)
	uim1 = mu_[0]; // left boundary condition
      if (i>0)
	uim1 = u[i-1];

      rhs[i] = ( 0.5*(uim1*uim1 - u[i]*u[i]) )/dx_;
      i++;
    }

    for (int i=0; i<NumMyElem_; ++i){
      rhs[i] += mu_[1]*exp(mu_[2] * (*xGrid_)[i]);
    }  
  }//end residual

  
  void jacobian(const state_type & u,
		jacobian_type & jac,
		const scalar_type /*t*/)
  {
    // //evaluate jacobian
    // if (jac.rows() == 0 || jac.cols()==0 ){
    //   jac.resize(u.size(), u.size());
    // }
    // // typedef Eigen::Triplet<double> Tr;
    // // std::vector<Tr> tripletList;
    // // tripletList.push_back( Tr(0,0,v_ij) );
    
    // jac = jacobian_type::Zero(jac.rows(), jac.cols());
    // jac(0,0) = -dxInv_ * u(0)*u(0);
    // for (int i=1; i<Ncell_; ++i){
    //   jac(i,i-1) = dxInv_ * u(i-1)*u(i-1);
    //   jac(i,i) = -dxInv_ * u(i-1)*u(i-1);     
    // }    
  }//end jacobian
  
private:  
  std::vector<scalar_type> mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain 
  const scalar_type xR_ = 100.0; // right side of domain
  int Ncell_; // # of cells
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size

  Epetra_MpiComm * comm_;
  rcp<Epetra_Map> dataMap_;
  int myRank_;
  int NumMyElem_;
  std::vector<int> myGel_;

  // mesh points coordinates
  rcp<nativeVec> xGrid_; 
  
  rcp<nativeVec> U_; // state vector
  rcp<nativeVec> U0_; // initial state vector
};

}//end namespace apps
#endif 
