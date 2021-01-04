/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_tpetra.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_DISCRETE_TIME_API_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_DISCRETE_TIME_API_HPP_

#include <numeric>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Tpetra_Core.hpp>

namespace pressio{ namespace apps{

class Burgers1dTpetraDiscreteTimeApi{
protected:
  using map_t		= Tpetra::Map<>;
  using nativeVec	= Tpetra::Vector<>;

  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  template<typename T> using stdrcp = std::shared_ptr<T>;
  using crs_graph_type = Tpetra::CrsGraph<>;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= nativeVec;
  using discrete_time_residual_type	= state_type;
  using jacobian_type	= Tpetra::CrsMatrix<>;
  using dense_matrix_type = Tpetra::MultiVector<>;

public:
  Burgers1dTpetraDiscreteTimeApi(std::vector<scalar_type> params,
		  int Ncell,
		  rcpcomm_t comm)
    : mu_{params}, Ncell_{Ncell}, comm_{comm}{
      this->setup();
    }

  Burgers1dTpetraDiscreteTimeApi() = delete;
  ~Burgers1dTpetraDiscreteTimeApi() = default;

public:
  rcpmap_t getDataMap(){
    return dataMap_;
  };

  state_type const & getInitialState() const{
    return *U0_;
  };

  discrete_time_residual_type createDiscreteTimeResidual() const{
    discrete_time_residual_type R(dataMap_);
    return R;
  }

  // computes: A = Jac B where B is dense
  dense_matrix_type createApplyDiscreteTimeJacobianResult(const dense_matrix_type & B) const
  {
    dense_matrix_type C( dataMap_, B.getNumVectors() );
    return C;
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
  			    Args && ... states) const
  {
    this->discreteTimeResidualImpl(step, time, dt, 
      R,std::forward<Args>(states)... );
  }


  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
  				 Args && ... states) const
  {
    this->applyDiscreteTimeJacobianImpl(step, time, dt, 
      B, A, std::forward<Args>(states)...);
  }

private:
  // implicit euler 
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    auto u_v = yn.getData();
    auto u_vnm1 = ynm1.getData();

    auto rhs_v = R.getDataNonConst();

    double valueFromLeft = 0.0;
    constexpr int tag_ = 1;
    if( myRank_ < totRanks_ - 1 ){
      MPI_Send( &u_v[NumMyElem_-1], 1, MPI_DOUBLE,
		myRank_+1, tag_, *comm_->getRawMpiComm() );
    }
    if( myRank_ > 0 ){
      MPI_Status status;
      MPI_Recv(&valueFromLeft, 1, MPI_DOUBLE,
	       myRank_-1, tag_,
	       *comm_->getRawMpiComm(), &status);
    }

    lo_t i=0;
    scalar_type uim1;
    for (auto const & it : myGel_){
      uim1 = valueFromLeft;
      if (it==0)
	uim1 = mu_[0]; // left boundary condition
      if (i>0)
	uim1 = u_v[i-1];

      rhs_v[i] = ( 0.5*(uim1*uim1 - u_v[i]*u_v[i]) )/dx_;
      i++;
    }

    auto xgrid_v = xGrid_->getDataNonConst();
    for (i=0; i<NumMyElem_; ++i){
      rhs_v[i] += mu_[1]*exp(mu_[2] * xgrid_v[i]);
      rhs_v[i] = u_v[i] - u_vnm1[i] - dt*rhs_v[i];
    }
  }

  // BDF2
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1,
                                const state_type & ynm2) const
  {
    auto u_v = yn.getData();
    auto u_vnm1 = ynm1.getData();
    auto u_vnm2 = ynm2.getData();

    auto rhs_v = R.getDataNonConst();

    double valueFromLeft = 0.0;
    constexpr int tag_ = 1;
    if( myRank_ < totRanks_ - 1 ){
      MPI_Send( &u_v[NumMyElem_-1], 1, MPI_DOUBLE,
		myRank_+1, tag_, *comm_->getRawMpiComm() );
    }
    if( myRank_ > 0 ){
      MPI_Status status;
      MPI_Recv(&valueFromLeft, 1, MPI_DOUBLE,
	       myRank_-1, tag_,
	       *comm_->getRawMpiComm(), &status);
    }

    lo_t i=0;
    scalar_type uim1;
    for (auto const & it : myGel_){
      uim1 = valueFromLeft;
      if (it==0)
	uim1 = mu_[0]; // left boundary condition
      if (i>0)
	uim1 = u_v[i-1];

      rhs_v[i] = ( 0.5*(uim1*uim1 - u_v[i]*u_v[i]) )/dx_;
      i++;
    }

    auto xgrid_v = xGrid_->getDataNonConst();
    for (i=0; i<NumMyElem_; ++i){
      rhs_v[i] += mu_[1]*exp(mu_[2] * xgrid_v[i]);
      rhs_v[i] = u_v[i] - 4./3.*u_vnm1[i] + 1./3.*u_vnm2[i] - 2./3.*dt*rhs_v[i];
    }
  }

  template <typename step_t, typename state_t>
  void applyDiscreteTimeJacobianImpl(const step_t & step,
				     const scalar_type & time,
				     const scalar_type & dt,
				     const dense_matrix_type & B,
				     dense_matrix_type & A,
				     const state_t & yn,
				     const state_t & ynm1) const

  {
    computeJacobianReplace(yn,ynm1, *Jac_,dt, time);
    Jac_->apply(B, A);
  }

  template <typename step_t, typename state_t>
  void applyDiscreteTimeJacobianImpl(const step_t & step,
				     const scalar_type & time,
				     const scalar_type & dt,
				     const dense_matrix_type & B,
				     dense_matrix_type & A,
				     const state_t & yn,
				     const state_t & ynm1,
                                     const state_t & ynm2) const

  {
    computeJacobianReplace(yn,ynm1,ynm2,*Jac_,dt, time);
    Jac_->apply(B, A);
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
    U0_ = std::make_shared<nativeVec>(dataMap_);
    U0_->putScalar(1.0);

    U_ = std::make_shared<nativeVec>(dataMap_);
    U_->putScalar(1.0);

    Teuchos::RCP<crs_graph_type> dataGraph =
      Teuchos::rcp(new crs_graph_type(dataMap_, nonZrPerRow_));

    this->assembleGraph(*dataGraph);
    Jac_ = std::make_shared<jacobian_type>(dataGraph);
    Jac_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);
  };

  void computeJacobianReplace(const state_type & u,
                              const state_type & unm1,
			      jacobian_type & jac,
                              const scalar_type dt,
			      const scalar_type /*t*/) const
  {
    const scalar_type buffer = exchangeFlux(u);
    auto u_v = u.getData();

    // resume fill because below we change entries
    jac.resumeFill();

    using tarr_dt = Teuchos::ArrayView<scalar_type>;
    using tarr_it = Teuchos::ArrayView<go_t>;
    std::array<scalar_type,1> va1;
    std::array<go_t,1> ci1;
    std::array<scalar_type,2> va2;
    std::array<go_t,2> ci2;

    for (lo_t i=0; i<NumMyElem_; i++){
      auto thisGID = myGel_[i];
      if (thisGID==0){
	va1[0] = dt*dxInv_ * u_v[0] + 1;
	ci1[0] = 0;
	jac.replaceGlobalValues(thisGID, tarr_it(ci1.data(),1),
				tarr_dt(va1.data(),1) );
      }
      else{
	ci2[0] = thisGID-1;
	ci2[1] = thisGID;

	if (i==0) va2[0] = -dt*dxInv_ * buffer;
	if (i>0) va2[0] = -dt*dxInv_ * u_v[i-1];

	va2[1] = -dt*dxInv_ * u_v[i] + 1;
	jac.replaceGlobalValues(thisGID, tarr_it(ci2.data(),2),
				tarr_dt(va2.data(),2) );
      }
    }

    jac.fillComplete();
  }//end

  void computeJacobianReplace(const state_type & u,
                              const state_type & unm1,
                              const state_type & unm2,
			      jacobian_type & jac,
                              const scalar_type dt,
			      const scalar_type /*t*/) const
  {
    const scalar_type buffer = exchangeFlux(u);
    auto u_v = u.getData();

    // resume fill because below we change entries
    jac.resumeFill();

    using tarr_dt = Teuchos::ArrayView<scalar_type>;
    using tarr_it = Teuchos::ArrayView<go_t>;
    std::array<scalar_type,1> va1;
    std::array<go_t,1> ci1;
    std::array<scalar_type,2> va2;
    std::array<go_t,2> ci2;

    for (lo_t i=0; i<NumMyElem_; i++){
      auto thisGID = myGel_[i];
      if (thisGID==0){
	va1[0] = 2./3.*dt*dxInv_ * u_v[0] + 1;
	ci1[0] = 0;
	jac.replaceGlobalValues(thisGID, tarr_it(ci1.data(),1),
				tarr_dt(va1.data(),1) );
      }
      else{
	ci2[0] = thisGID-1;
	ci2[1] = thisGID;

	if (i==0) va2[0] = -2./3.*dt*dxInv_ * buffer;
	if (i>0) va2[0] = -2./3.*dt*dxInv_ * u_v[i-1];

	va2[1] = 2./3.*dt*dxInv_ * u_v[i] + 1;
	jac.replaceGlobalValues(thisGID, tarr_it(ci2.data(),2),
				tarr_dt(va2.data(),2) );
      }
    }

    jac.fillComplete();
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
    auto u_v = u.getData();

    // to populate the jacobian each process needs the last grid
    // point solution from the previous process
    scalar_type buffer {0.0};
    MPI_Status st;
    constexpr int tag_ = 1;
    if (myRank_ < totRanks_-1){
      MPI_Send(&u_v[NumMyElem_-1], 1, MPI_DOUBLE,
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

  mutable stdrcp<nativeVec> U_{}; // state vector
  mutable stdrcp<nativeVec> U0_{}; // initial state vector
  stdrcp<jacobian_type> Jac_{};

};//end class

}} //namespace pressio::apps
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_TPETRA_HPP_
