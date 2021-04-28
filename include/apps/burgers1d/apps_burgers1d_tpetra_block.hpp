/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_tpetra_block.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_BLOCK_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_BLOCK_HPP_

#include <numeric>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockVector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Tpetra_Core.hpp>

namespace pressio{ namespace apps{

class Burgers1dTpetraBlock{
protected:
  using map_t		= Tpetra::Map<>;
  using nativeBlockVec	= Tpetra::BlockVector<>;
  using nativeVec	= Tpetra::Vector<>;
  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  template<typename T> using stdrcp = std::shared_ptr<T>;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= nativeBlockVec;
  using velocity_type	= state_type;
  using jacobian_type = Tpetra::BlockCrsMatrix<>;
  using crs_graph_type = Tpetra::CrsGraph<>;

public:
  Burgers1dTpetraBlock(std::vector<scalar_type> params,
		       int Ncell,
		       rcpcomm_t comm)
    : mu_{params}, Ncell_{Ncell}, comm_{comm}{
      this->setup();
    }

public:
  rcpmap_t getDataMap(){
    return dataMap_;
  };

  state_type const & getInitialState() const{
    return *U0_;
  };

  void velocity(const state_type & u,
		const scalar_type t,
		velocity_type & rhs) const
  {

    double valueFromLeft = 0.0;
    constexpr int tag_ = 1;

    if( myRank_ < totRanks_ - 1 ){
      auto uC = u.getLocalBlock(NumMyElem_ - 1);
      MPI_Send( &uC(0), 1, MPI_DOUBLE,
		myRank_+1, tag_, *comm_->getRawMpiComm() );
    }
    if( myRank_ > 0 ){
      MPI_Status status;
      MPI_Recv(&valueFromLeft, 1, MPI_DOUBLE,
	       myRank_-1, tag_,
	       *comm_->getRawMpiComm(), &status);    }
    lo_t i=0;
    scalar_type uim1;
    scalar_type rhsLoc;
    for (auto const & it : myGel_)
    {
      uim1 = valueFromLeft;
      const auto uCi = u.getLocalBlock(i)(0);
      if (it==0)
	uim1 = mu_[0];
      if (i>0)
	uim1 = u.getLocalBlock(i-1)(0);

      rhsLoc = ( 0.5*(uim1*uim1 - uCi*uCi) )/dx_;
      rhs.replaceLocalValues(i,&rhsLoc);
      i++;
    }

    auto xgrid_v = xGrid_->getDataNonConst();
    for (i=0; i<NumMyElem_; ++i){
      auto rhsVal = rhs.getLocalBlock(i);
      rhsLoc = rhsVal(0) + mu_[1]*exp(mu_[2] * xgrid_v[i]);
      rhs.replaceLocalValues(i,&rhsLoc);
    }
  }

  // computes: A = Jac B where B is dense
  void applyJacobian(const state_type & y,
         const Tpetra::BlockMultiVector<> & B,
         scalar_type t,
         Tpetra::BlockMultiVector<> & A) const
  {
    computeJacobianUpdate(y, *Jac_, t);
    const auto B_vv = B.getMultiVectorView();
    auto A_vv     = A.getMultiVectorView();
    Jac_->apply(B_vv, A_vv);
    //Jac_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);
  }

  velocity_type createVelocity() const
  {
    velocity_type R(*dataMap_,1);
    return R;
  }

  // computes: A = Jac B where B is dense
  Tpetra::BlockMultiVector<>
  createApplyJacobianResult(const Tpetra::BlockMultiVector<> & B) const
  {
    Tpetra::BlockMultiVector<> A( *dataMap_, blockSize, B.getNumVectors() );
    return A;
  }

protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();
    // distribute cells
    dataMap_ = Teuchos::rcp(new map_t(Ncell_, 0, comm_));

    NumMyElem_ = dataMap_->getNodeNumElements();
    auto minGId = dataMap_->getMinGlobalIndex();

    myGel_.resize(NumMyElem_);
    std::iota(myGel_.begin(), myGel_.end(), minGId);
    //dataMap_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    if (myRank_==1){
      for (auto it : myGel_)
	std::cout << it << std::endl;
    }

    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;

    // grid
    xGrid_ = std::make_shared<nativeVec>(dataMap_);
    auto xGridv = xGrid_->getDataNonConst();
    lo_t i=0;
    for (auto const & it : myGel_){
      xGridv[i] = dx_*it + dx_*0.5;
      i++;
    };
    //xGrid_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    // init condition
    U0_ = std::make_shared<nativeBlockVec>(*dataMap_,blockSize);
    U0_->putScalar(1.0);

    U_ = std::make_shared<nativeBlockVec>(*dataMap_,blockSize);
    U_->putScalar(1.0);

    // construct a graph for the block matrix
    crs_graph_type dataGraph(dataMap_,nonZrPerRow_);
    this->assembleGraph(dataGraph);
    Jac_ = std::make_shared<jacobian_type>(dataGraph,blockSize);

    //this->computeJacobianUpdate(*U_, *Jac_, 0.0);
    Jac_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);
  };

  void computeJacobianUpdate(const state_type & u,
			     jacobian_type & jac,
			     const scalar_type ) const
  {
    const scalar_type buffer = exchangeFlux(u);

    const lo_t* localColInds;
    scalar_type* vals;
    for (lo_t i=0; i<NumMyElem_; i++)
    {
      auto thisGID = myGel_[i];
      auto uC = u.getLocalBlock(i);
      if (thisGID==0)
	{
	  lo_t numEntries = 1;
	  int err = 0;
	  err = jac.getLocalRowView(i, localColInds, vals, numEntries);
	  if (err != 0) {
	    break;
	  }
	  vals[0]  =  -dxInv_ * uC(0);
	}
      else
	{
	  lo_t numEntries = 2;
	  int err = 0;
	  err = jac.getLocalRowView(i, localColInds, vals, numEntries);
    (void)err;
    
	  if (i==0){
	    // don't know why the indices here must be inverted, but this is right/
	    // looks like block does strange here for storing the col indices of this point
	    vals[1] = dxInv_ * buffer;
	    vals[0] = -dxInv_ * uC(0);
	  }
	  if (i>0){
	    vals[0] = dxInv_ * u.getLocalBlock(i-1)(0);
	    vals[1] = -dxInv_ * uC(0);
	  }
	}
    }
  }//end

  void assembleGraph(crs_graph_type & graph)
  {
    // using tarr_dt = Teuchos::ArrayView<scalar_type>;
    using tarr_it = Teuchos::ArrayView<go_t>;
    std::array<go_t,1> ci1;
    std::array<go_t,2> ci2;
    for (lo_t i=0; i<NumMyElem_; i++){
      auto thisGID = myGel_[i];
      if (thisGID==0){
	ci1[0] = 0;
	graph.insertGlobalIndices(thisGID, tarr_it(ci1.data(),1));
      }
      else{
	ci2[0] = thisGID-1;
	ci2[1] = thisGID;
	graph.insertGlobalIndices(thisGID, tarr_it(ci2.data(),2));
      }
    }
    graph.fillComplete();
  }//end

  scalar_type exchangeFlux(const state_type & u) const
  {
    // to populate the jacobian each process needs the last grid
    // point solution from the previous process
    scalar_type buffer {0.0};
    MPI_Status st;
    constexpr int tag_ = 1;
    if (myRank_ < totRanks_-1){
      auto uC = u.getLocalBlock(NumMyElem_ - 1);
      MPI_Send(&uC(0), 1, MPI_DOUBLE,
	       myRank_+1, tag_, *comm_->getRawMpiComm());
    }
    if (myRank_ >= 1){
      MPI_Recv(&buffer, 1, MPI_DOUBLE, myRank_-1, tag_,
	       *comm_->getRawMpiComm(), &st);
    }
    return buffer;
  }

protected:
  std::vector<scalar_type> mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain
  const scalar_type xR_ = 100.0; // right side of domain
  const lo_t blockSize = 1;
  int Ncell_{}; // # of cells
  scalar_type dx_{}; // cell size
  scalar_type dxInv_{}; // inv of cell size
  const int nonZrPerRow_ = 2;
  stdrcp<nativeVec> xGrid_{}; // mesh points coordinates

  Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
  rcpcomm_t comm_{};
  rcpmap_t dataMap_{};

  int myRank_{};
  int totRanks_{};
  lo_t NumMyElem_{};
  std::vector<go_t> myGel_{};

  mutable stdrcp<nativeBlockVec> U_{}; // state vector
  mutable stdrcp<nativeBlockVec> U0_{}; // initial state vector
  stdrcp<jacobian_type> Jac_{};

};//end class

}} //namespace pressio::apps
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_BLOCK_HPP_
