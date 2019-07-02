
#include "apps_steady_linear_adv_diff_2d_epetra.hpp"

#ifdef HAVE_TRILINOS
namespace pressio{ namespace apps{

void SteadyLinAdvDiff2dEpetra::createMap(){
  // total number of dof (we only consider the interior points)
  numGlobalDof_ = NxDof_ * NyDof_;
  myMap_ = std::make_shared<Epetra_Map>(numGlobalDof_, 0, comm_);
}

void SteadyLinAdvDiff2dEpetra::setup(){
  createMap();
  dofPerProc_= myMap_->NumMyElements();
  MyGlobalDof_ = myMap_->MyGlobalElements();

  x_ = std::make_shared<nativeVec>(*myMap_);
  y_ = std::make_shared<nativeVec>(*myMap_);
  u_ = std::make_shared<nativeVec>(*myMap_);
  v_ = std::make_shared<nativeVec>(*myMap_);
  int gi{}; int gj{};
  for (auto i = 0; i<dofPerProc_; i++){
    auto GID = MyGlobalDof_[i];
    globalIDToGiGj(GID, gi, gj);
    (*x_)[i] = dx_ + gi * dx_;
    (*y_)[i] = gj * dy_;
    auto xval = (*x_)[i];
    auto yval = (*y_)[i];
    (*u_)[i] = -std::sin(M_PI*xval) * std::cos(M_PI*yval);
    (*v_)[i] = std::cos(M_PI*xval) * std::sin(M_PI*yval);
  }

  A_ = std::make_shared<nativeMatrix>(Copy, *myMap_, maxNonZeroPerRow_);
  T_ = std::make_shared<nativeVec>(*myMap_);
  f_ = std::make_shared<nativeVec>(*myMap_);
}


void SteadyLinAdvDiff2dEpetra::assembleMatrix() const{

  int gi{}; int gj{};
  int k{};
  std::array<int, 5> colInd{};
  std::array<scalar_type, 5> values{};
  scalar_type c_T_ip1{};
  scalar_type c_T_im1{};
  scalar_type c_T_ij {};
  scalar_type c_T_jp1{};
  scalar_type c_T_jm1{};
  int gid_ip1{};
  int gid_im1{};
  int gid_jp1{};
  int gid_jm1{};

  for (auto i=0; i<dofPerProc_; i++)
  {
    // global ID of current point
    auto GID = MyGlobalDof_[i];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    const scalar_type uij = (*u_)[i];
    const scalar_type vij = (*v_)[i];

    k = 0;
    c_T_ij    = - invPrRe_*dxSqInv_ - invPrRe_*dySqInv_;
    values[k] = 2.0*c_T_ij;
    colInd[k] = GID;
    k++;

    if (gi>=1){
      c_T_im1 = invPrRe_*dxSqInv_ + uij*dx2Inv_;
      gid_im1 = GID - 1;
      values[k] = c_T_im1;
      colInd[k] = gid_im1;
      k++;
    }

    if (gi<NxDof_-1){
      c_T_ip1 = invPrRe_*dxSqInv_ - uij*dx2Inv_;
      gid_ip1 = GID + 1;
      values[k] = c_T_ip1;
      colInd[k] = gid_ip1;
      k++;
    }

    if (gj>=1 and gj<NyDof_-1){
      c_T_jm1 = invPrRe_*dySqInv_ + vij*dy2Inv_;
      gid_jm1 = GID - NxDof_;
      values[k] = c_T_jm1;
      colInd[k] = gid_jm1;
      k++;

      c_T_jp1 = invPrRe_*dySqInv_ - vij*dy2Inv_;
      gid_jp1 = GID + NxDof_;
      values[k] = c_T_jp1;
      colInd[k] = gid_jp1;
      k++;
    }

    if (gj==0){
      c_T_jp1 = 2.0*invPrRe_*dySqInv_;
      gid_jp1 = GID + NxDof_;
      values[k] = c_T_jp1;
      colInd[k] = gid_jp1;
      k++;
    }

    if (gj==NyDof_-1){
      c_T_jm1 = 2.0*invPrRe_*dySqInv_;
      gid_jm1 = GID - NxDof_;
      values[k] = c_T_jm1;
      colInd[k] = gid_jm1;
      k++;
    }

    if (A_->Filled())
      A_->ReplaceGlobalValues(GID, k,
			      values.data(), colInd.data());
    else
      A_->InsertGlobalValues(GID, k,
			     values.data(), colInd.data());
  }//loop

  if(!A_->Filled())
    A_->FillComplete();
}


void SteadyLinAdvDiff2dEpetra::fillRhs() const{
  f_->PutScalar( static_cast<scalar_type>(0) );

  int gi{}; int gj{};
  for (auto i=0; i<dofPerProc_; i++)
  {
    // global ID of current point
    auto GID = MyGlobalDof_[i];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    if (gi==0){
      auto value = invPrRe_*dxSqInv_ + (*u_)[i]*dx2Inv_;
      (*f_)[i] = -value * Tleft_;
    }

    if (gi==NxDof_-1){
      auto value = invPrRe_*dxSqInv_ - (*u_)[i]*dx2Inv_;
      (*f_)[i] = -value * Tright_;
    }
  }//loop
}


void SteadyLinAdvDiff2dEpetra::solve(){
  Epetra_LinearProblem Problem(A_.get(), T_.get(), f_.get());
  AztecOO Solver(Problem);
  Solver.Iterate(500, solveTolerance_);
  Solver.NumIters();
  Solver.TrueResidual();
}

}} //namespace pressio::apps
#endif
