
#if HAVE_TRILINOS
#include "apps_burgers1d_tpetra.hpp"

namespace rompp{ namespace apps{

void Burgers1dTpetra::residual(const state_type & u,
			       residual_type & rhs,
			       const scalar_type /* t */) const
{
  auto u_v = u.getData();
  auto rhs_v = rhs.getDataNonConst();

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
  }
}
//-------------------------------------------------------


Burgers1dTpetra::scalar_type
Burgers1dTpetra::exchangeFlux(const state_type & u) const{

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
//-------------------------------------------------------


void Burgers1dTpetra::computeJacobianInsert(const state_type & u,
					    jacobian_type & jac,
					    const scalar_type /*t*/) const{

  const scalar_type buffer = exchangeFlux(u);
  auto u_v = u.getData();

  using tarr_dt = Teuchos::ArrayView<scalar_type>;
  using tarr_it = Teuchos::ArrayView<go_t>;
  std::array<scalar_type,1> va1;
  std::array<go_t,1> ci1;
  std::array<scalar_type,2> va2;
  std::array<go_t,2> ci2;

  for (lo_t i=0; i<NumMyElem_; i++){
    auto thisGID = myGel_[i];
    if (thisGID==0){
	va1[0] = -dxInv_ * u_v[0];
	ci1[0] = 0;
	jac.insertGlobalValues(thisGID, tarr_it(ci1.data(),1),
			       tarr_dt(va1.data(),1) );
      }
    else{
      ci2[0] = thisGID-1;
      ci2[1] = thisGID;
      if (i==0) va2[0] = dxInv_ * buffer;
      if (i>0) va2[0] = dxInv_ * u_v[i-1];

      va2[1] = -dxInv_ * u_v[i];
      jac.insertGlobalValues(thisGID, tarr_it(ci2.data(),2),
			     tarr_dt(va2.data(),2) );
    }
  }

  jac.fillComplete();
}//end jacobian



void Burgers1dTpetra::computeJacobianReplace(const state_type & u,
					    jacobian_type & jac,
					    const scalar_type /*t*/) const{
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
	va1[0] = -dxInv_ * u_v[0];
	ci1[0] = 0;
	jac.replaceGlobalValues(thisGID, tarr_it(ci1.data(),1),
				tarr_dt(va1.data(),1) );
      }
    else{
      ci2[0] = thisGID-1;
      ci2[1] = thisGID;
      if (i==0) va2[0] = dxInv_ * buffer;
      if (i>0) va2[0] = dxInv_ * u_v[i-1];

      va2[1] = -dxInv_ * u_v[i];
      jac.replaceGlobalValues(thisGID, tarr_it(ci2.data(),2),
			      tarr_dt(va2.data(),2) );
    }
  }

  jac.fillComplete();
}//end jacobian replace


}} //namespace rompp::apps
#endif
