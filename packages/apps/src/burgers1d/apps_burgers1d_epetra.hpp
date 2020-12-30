/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_epetra.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_HPP_

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

namespace pressio{ namespace apps{

class Burgers1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;

/* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using velocity_type	= state_type;
  using jacobian_type	= Epetra_CrsMatrix;

public:
  Burgers1dEpetra(std::vector<scalar_type> params,
		  int Ncell, Epetra_MpiComm * comm)
    : mu_(params), Ncell_(Ncell), comm_(comm){
    this->setup();
  }

  Burgers1dEpetra() = delete;
  ~Burgers1dEpetra() = default;

public:
  Epetra_Map const & getDataMap(){
    return *dataMap_;
  };

  state_type const & getInitialState() const{
    return *U0_;
  };

  velocity_type createVelocity() const{
    Epetra_Vector R(*dataMap_);
    return R;
  }

  // computes: A = Jac B where B is a multivector
  Epetra_MultiVector createApplyJacobianResult(const Epetra_MultiVector & B) const{
    Epetra_MultiVector C( Jac_->RangeMap(), B.NumVectors() );
    return C;
  }

  void velocity(const state_type & u,
		const scalar_type /* t */,
		velocity_type & rhs) const
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
  }//end velocity
  //-------------------------------------------------------

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const Epetra_MultiVector & B,
		     scalar_type t,
		     Epetra_MultiVector & A) const{
    assert( Jac_->NumGlobalCols() == B.GlobalLength() );
    assert( A.GlobalLength() == Jac_->NumGlobalRows() );
    assert( A.NumVectors() == B.NumVectors() );
    // compute jacobian
    jacobian(y, t, *Jac_);
    Jac_->Print(std::cout);
    // multiply
    Jac_->Multiply(false, B, A);
    //A.Print(std::cout);
  }


  void jacobian(const state_type & u,
		const scalar_type /*t*/,
		jacobian_type & jac) const
  {
    // to populate the jacobian each process needs the last grid
    // point solution from the previous process
    double buffin {0.0};
    MPI_Status st{};
    MPI_Comm mpiComm = comm_->Comm();
    if (myRank_ < totRanks_-1){
      double tosend = u[NumMyElem_-1];
      MPI_Send(&tosend, 1, MPI_DOUBLE, myRank_+1, 1, mpiComm);
    }
    if (myRank_ >= 1){
      MPI_Recv(&buffin, 1, MPI_DOUBLE, myRank_-1, 1, mpiComm, &st);
    }

    std::vector<int> Indices {0, 0};
    std::vector<double> Values {0., 0.};
    for (int i=0; i<NumMyElem_; i++)
      {
	int thisGID = myGel_[i]; // global ID
	if (thisGID==0){
	  Indices[0] = 0;
	  Values[0] = -dxInv_ * u[0];
	  if (jac.Filled())
	    jac.ReplaceGlobalValues(thisGID, 1,
				    Values.data(), Indices.data() );
	  else
	    jac.InsertGlobalValues(thisGID, 1,
				   Values.data(), Indices.data() );
	}
	else{
	  Indices[0] = thisGID-1;
	  Indices[1] = thisGID;
	  if (i==0){
	    Values[0] = dxInv_ * buffin;
	    Values[1] = -dxInv_ * u[i];
	  }
	  if (i>0){
	    Values[0] = dxInv_ * u[i-1];
	    Values[1] = -dxInv_ * u[i];
	  }

	  if (jac.Filled())
	    jac.ReplaceGlobalValues(thisGID, 2,
				    Values.data(), Indices.data() );
	  else
	    jac.InsertGlobalValues(thisGID, 2,
				   Values.data(), Indices.data() );
	}
      }
    if (!jac.Filled())
      jac.FillComplete();
  }//end jacobian
  //-------------------------------------------------------

protected:
  void setup(){
    // distribute cells
    dataMap_ = std::make_shared<Epetra_Map>(Ncell_,0,*comm_);
    dataMap_->Print(std::cout);

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
    Jac_ = std::make_shared<Epetra_CrsMatrix>(Copy, *dataMap_, nonZrPerRow_);
  };

protected:
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

  mutable rcp<nativeVec> U_; // state vector
  mutable rcp<nativeVec> U0_; // initial state vector
  std::shared_ptr<Epetra_CrsMatrix> Jac_;
};//end class

}} //namespace pressio::apps

#endif  // APPS_BURGERS1D_APPS_BURGERS1D_EPETRA_HPP_
