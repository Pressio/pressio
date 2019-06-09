
#if HAVE_TRILINOS
#include "apps_burgers1d_epetra.hpp"

namespace rompp{ namespace apps{

void Burgers1dEpetra::residual(const state_type & u,
			       residual_type & rhs,
			       const scalar_type /* t */) const
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
//-------------------------------------------------------

void Burgers1dEpetra::jacobian(const state_type & u,
			       jacobian_type & jac,
			       const scalar_type /*t*/) const{

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

}} //namespace rompp::apps
#endif
