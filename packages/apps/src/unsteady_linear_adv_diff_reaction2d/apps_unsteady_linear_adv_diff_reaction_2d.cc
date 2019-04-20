
#include "apps_unsteady_linear_adv_diff_reaction_2d.hpp"

namespace rompp{ namespace apps{

void UnsteadyLinAdvDiffReac2dEpetra::createGridMap(){
  // total number of unknown grid points (we only consider the interior points)
  numGlobalUnkGpt_ = Nx_ * Ny_;
  gridMap_ = std::make_shared<Epetra_Map>(numGlobalUnkGpt_, 0, comm_);
}

void UnsteadyLinAdvDiffReac2dEpetra::setupPhysicalGrid(){
  gptPerProc_	  = gridMap_->NumMyElements();
  MyGlobalUnkGpt_ = gridMap_->MyGlobalElements();

  x_ = std::make_shared<nativeVec>(*gridMap_);
  y_ = std::make_shared<nativeVec>(*gridMap_);
  u_ = std::make_shared<nativeVec>(*gridMap_);
  v_ = std::make_shared<nativeVec>(*gridMap_);
  int gi{}; int gj{};
  for (auto i = 0; i<gptPerProc_; i++){
    auto GID = MyGlobalUnkGpt_[i];
    //std::cout << i << " " << GID << std::endl;
    globalIDToGiGj(GID, gi, gj);
    (*x_)[i] = dx_ + gi * dx_;
    (*y_)[i] = gj * dy_;
    auto xval = (*x_)[i];
    auto yval = (*y_)[i];
    (*u_)[i] = -std::sin(M_PI*xval) * std::cos(M_PI*yval);
    (*v_)[i] = std::cos(M_PI*xval) * std::sin(M_PI*yval);
  }
}

void UnsteadyLinAdvDiffReac2dEpetra::createSolutionMap(){
  // total number of dof (we only consider the interior points)
  // we need to account for the fact that we have 3 fields per grid
  numGlobalDof_ = this_t::numSpecies_ * Nx_ * Ny_;
  dofMap_ = std::make_shared<Epetra_Map>(numGlobalDof_, 0, comm_);
}

void UnsteadyLinAdvDiffReac2dEpetra::setupFields(){
  constexpr auto maxNonZ = this_t::maxNonZeroPerRow_;

  dofPerProc_= dofMap_->NumMyElements();
  MyGlobalDof_ = dofMap_->MyGlobalElements();

  s1_ = std::make_shared<nativeVec>(*gridMap_);
  s2_ = std::make_shared<nativeVec>(*gridMap_);
  s3_ = std::make_shared<nativeVec>(*gridMap_);

  A_ = std::make_shared<nativeMatrix>(Copy, *dofMap_, maxNonZ);
  chemForcing_ = std::make_shared<nativeVec>(*dofMap_);
  state_ = std::make_shared<nativeVec>(*dofMap_);
  A_->PutScalar(0);
  chemForcing_->PutScalar(0);
  state_->PutScalar(0);
}


void UnsteadyLinAdvDiffReac2dEpetra::fillSource1(){
  for (auto i=0; i<gptPerProc_; i++){
    const scalar_type xij = (*x_)[i];
    const scalar_type yij = (*y_)[i];
    const scalar_type dx = (xij - oPtS1[0]);
    const scalar_type dy = (yij - oPtS1[1]);
    const scalar_type distance = dx*dx + dy*dy;
    (*s1_)[i] = ( std::sqrt(distance) <= rS1) ? 0.1 : 0.0;
  }
}

void UnsteadyLinAdvDiffReac2dEpetra::fillSource2(){
  for (auto i=0; i<gptPerProc_; i++){
    const scalar_type xij = (*x_)[i];
    const scalar_type yij = (*y_)[i];
    const scalar_type dx = (xij - oPtS2[0]);
    const scalar_type dy = (yij - oPtS2[1]);
    const scalar_type distance = dx*dx + dy*dy;
    (*s2_)[i] = ( std::sqrt(distance) <= rS2) ? 0.1 : 0.0;
  }
}

void UnsteadyLinAdvDiffReac2dEpetra::fillSource3(){
  s3_->PutScalar(0);
}


void UnsteadyLinAdvDiffReac2dEpetra::setup(){
  createGridMap();
  createSolutionMap();
  setupPhysicalGrid();
  setupFields();
  fillSource1();
  fillSource2();
}

void UnsteadyLinAdvDiffReac2dEpetra::assembleMatrix() const
{
  /* note that R is the residual vector, which has
   * a dofMap_ because it contains all dofs.
   * Once again, the dofs are the numFields * numOfUnkownGridPoints

   * Note also that this matrix computes convective and diffusion terms
   * assuming they are on the RHS
   So: df/dt = -conv + diffusion + else

   This assembles the FD matrix for the -conv + diffusion
   */

  constexpr auto maxNonZ = UnsteadyLinAdvDiffReac2dEpetra::maxNonZeroPerRow_;
  int gi{}; int gj{};
  int k{};
  std::array<int, maxNonZ> colInd{};
  std::array<scalar_type, maxNonZ> values{};

  int l = {0};

  // loop over grid points
  for (auto i=0; i<gptPerProc_; i++)
  {
    // global ID of current grid point
    auto GID = MyGlobalUnkGpt_[i];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);
    // the local valocity
    auto uij = (*u_)[i];
    auto vij = (*v_)[i];

    // for a given grid point, loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      values = {};
      colInd = {};

      // global ID of current UNKNOWN
      auto dofGID = MyGlobalDof_[l];

      // get the diffusivity for current species
      auto D = eps_;

      // add diagonal
      k = 0;
      values[k] = 2.0*(-D * dxSqInv_ - D * dySqInv_);
      colInd[k] = dofGID;
      k++;

      // i-1, j
      if (gi>=1){
      	values[k] = D*dxSqInv_ + uij*dx2Inv_;
      	colInd[k] = dofGID - numSpecies_;
      	k++;
      }

      // i+1, j
      if (gi<Nx_-1){
      	values[k] = D*dxSqInv_ - uij*dx2Inv_;
      	colInd[k] = dofGID + numSpecies_;
      	k++;
      }

      // i, j-1 and i, j+1
      if (gj>=1 and gj<Ny_-1){
      	values[k] = D*dySqInv_ + vij*dy2Inv_;
      	colInd[k] = dofGID - Nx_*numSpecies_;
      	k++;

      	values[k] = D*dySqInv_ - vij*dy2Inv_;
      	colInd[k] = dofGID + Nx_*numSpecies_;
      	k++;
      }

      // bottom wall we have homog Neumann BC
      if (gj==0){
      	values[k] = 2.0*D*dySqInv_;
      	colInd[k] = dofGID + Nx_*numSpecies_;
      	k++;
      }

      // top wall we have homog Neumann BC
      if (gj==Ny_-1){
      	values[k] = 2.0*D*dySqInv_;
      	colInd[k] = dofGID - Nx_*numSpecies_;
      	k++;
      }

      if (A_->Filled())
	A_->ReplaceGlobalValues(dofGID, k,
				values.data(), colInd.data());
      else
	A_->InsertGlobalValues(dofGID, k,
			       values.data(), colInd.data());

      // update counter
      l++;
    }//over dof
  }//loop over grid pts

  if(!A_->Filled())
    A_->FillComplete();
}// end method


void UnsteadyLinAdvDiffReac2dEpetra::computeChem
( const state_type & C ) const{

  chemForcing_->PutScalar( static_cast<scalar_type>(0) );
  unsigned int k = 0;
  // loop over grid points
  for (auto i=0; i<gptPerProc_; i++){
    // for a given grid point, loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++){
      if (iDof == 0)
  	(*chemForcing_)[k] = -K_*C[k] * C[k+1] + (*s1_)[i];
      else if (iDof == 1)
  	(*chemForcing_)[k] = -K_ * C[k-1] * C[k] + (*s2_)[i];
      else if (iDof == 2)
  	(*chemForcing_)[k] = K_*C[k-2]*C[k-1] - K_ * C[k] + (*s3_)[i];

      // update indexer
      k++;
    }
  }
}

}} //namespace rompp::apps
