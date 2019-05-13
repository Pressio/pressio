
#include "apps_unsteady_nonlinear_adv_diff_reaction_2d_epetra.hpp"

namespace rompp{ namespace apps{

void UnsteadyNonLinAdvDiffReac2dEpetra::createGridMap(){
  // total number of unknown grid points (we only consider the interior points)
  numGlobalUnkGpt_ = Nx_ * Ny_;
  gridMap_ = std::make_shared<Epetra_Map>(numGlobalUnkGpt_, 0, comm_);
  gptPerProc_	  = gridMap_->NumMyElements();
  MyGlobalUnkGpt_ = gridMap_->MyGlobalElements();
  //gridMap_->Print(std::cout);
}

void UnsteadyNonLinAdvDiffReac2dEpetra::setupPhysicalGrid(){
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

void UnsteadyNonLinAdvDiffReac2dEpetra::createSolutionMap(){
  // total number of dof (we only consider the interior points)
  // we need to account for the fact that we have 3 fields per grid
  // also, we do not want to break dofs, so if a rank owns
  // a specific grid point, then it owns all dofs associated
  // with that grid point

  // get the minimum GID of grid points I own
  const auto myMinGid = gridMap_->MinMyGID();

  // find the min GID of the dof map
  auto myMinDofGID = myMinGid * this_t::numSpecies_;

  numGlobalDof_ = this_t::numSpecies_ * Nx_ * Ny_;
  const auto dofPerProc_ = gptPerProc_*this_t::numSpecies_;

  std::vector<int> myGDofs(dofPerProc_);
  std::iota( myGDofs.begin(), myGDofs.end(), myMinDofGID);
  dofMap_ = std::make_shared<Epetra_Map>
    (numGlobalDof_, dofPerProc_, myGDofs.data(), 0, comm_);
  MyGlobalDof_ = dofMap_->MyGlobalElements();
  //dofMap_->Print(std::cout);
}

void UnsteadyNonLinAdvDiffReac2dEpetra::setupFields(){
  constexpr auto maxNonZ = this_t::maxNonZeroPerRow_;

  s1_ = std::make_shared<nativeVec>(*gridMap_);
  s2_ = std::make_shared<nativeVec>(*gridMap_);
  s3_ = std::make_shared<nativeVec>(*gridMap_);

  FDMat_ = std::make_shared<nativeMatrix>(Copy, *dofMap_, maxNonZ);
  chemForcing_ = std::make_shared<nativeVec>(*dofMap_);
  state_ = std::make_shared<nativeVec>(*dofMap_);
  FDMat_->PutScalar(0);
  chemForcing_->PutScalar(0);
  state_->PutScalar(0);
}

void UnsteadyNonLinAdvDiffReac2dEpetra::fillSource1(){
  for (auto i=0; i<gptPerProc_; i++){
    const scalar_type xij = (*x_)[i];
    const scalar_type yij = (*y_)[i];
    const scalar_type dx = (xij - oPtS1[0]);
    const scalar_type dy = (yij - oPtS1[1]);
    const scalar_type distance = dx*dx + dy*dy;
    (*s1_)[i] = ( std::sqrt(distance) <= rS1) ? 0.1 : 0.0;
  }
}

void UnsteadyNonLinAdvDiffReac2dEpetra::fillSource2(){
  for (auto i=0; i<gptPerProc_; i++){
    const scalar_type xij = (*x_)[i];
    const scalar_type yij = (*y_)[i];
    const scalar_type dx = (xij - oPtS2[0]);
    const scalar_type dy = (yij - oPtS2[1]);
    const scalar_type distance = dx*dx + dy*dy;
    (*s2_)[i] = ( std::sqrt(distance) <= rS2) ? 0.1 : 0.0;
  }
}

void UnsteadyNonLinAdvDiffReac2dEpetra::fillSource3(){
  s3_->PutScalar(0);
}

void UnsteadyNonLinAdvDiffReac2dEpetra::setup(){
  createGridMap();
  createSolutionMap();
  setupPhysicalGrid();
  setupFields();
  fillSource1();
  fillSource2();
}

void UnsteadyNonLinAdvDiffReac2dEpetra::assembleFDMatrix() const
{
  /* note that R is the residual vector, which has
   * a dofMap_ because it contains all dofs.
   * Once again, the dofs are the numFields * numOfUnkownGridPoints

   * Note also that this matrix computes convective and diffusion terms
   * assuming they are on the RHS
   So: df/dt = -conv + diffusion + else

   This assembles the FD matrix for the -conv + diffusion
   */

  constexpr auto maxNonZ = UnsteadyNonLinAdvDiffReac2dEpetra::maxNonZeroPerRow_;
  int gi{}; int gj{}; int k{};
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
      // global ID of current UNKNOWN
      auto dofGID = MyGlobalDof_[l];

      // add diagonal
      k = 0;
      values[k] = -2.0*(eps_*dxSqInv_ + eps_*dySqInv_);
      colInd[k] = dofGID;
      k++;

      // i-1, j
      if (gi>=1){
      	values[k] = eps_*dxSqInv_ + uij*dx2Inv_;
      	colInd[k] = dofGID - numSpecies_;
      	k++;
      }

      // i+1, j
      if (gi<Nx_-1){
      	values[k] = eps_*dxSqInv_ - uij*dx2Inv_;
      	colInd[k] = dofGID + numSpecies_;
      	k++;
      }

      // i, j-1 and i, j+1
      if (gj>=1 and gj<Ny_-1){
      	values[k] = eps_*dySqInv_ + vij*dy2Inv_;
      	colInd[k] = dofGID - Nx_*numSpecies_;
      	k++;

      	values[k] = eps_*dySqInv_ - vij*dy2Inv_;
      	colInd[k] = dofGID + Nx_*numSpecies_;
      	k++;
      }

      // bottom wall we have homog Neumann BC
      if (gj==0){
      	values[k] = 2.0*eps_*dySqInv_;
      	colInd[k] = dofGID + Nx_*numSpecies_;
      	k++;
      }

      // top wall we have homog Neumann BC
      if (gj==Ny_-1){
      	values[k] = 2.0*eps_*dySqInv_;
      	colInd[k] = dofGID - Nx_*numSpecies_;
      	k++;
      }

      // the following 3 terms are the ones for chemistry but
      // are set to zero here because are used later
      // to compute jacobian, they are not used for FD matrix.
      // However, we need to insert them here so that we can
      // replace their values in the matrix later
      // without epetra complaining.
      if (iDof == 0){
      	values[k] = 0.0;
      	colInd[k] = dofGID+1;
      	k++;
      }

      if (iDof == 1){
      	values[k] = 0.0;
      	colInd[k] = dofGID-1;
      	k++;
      }

      if (iDof == 2){
      	values[k] = 0.0;
	colInd[k] = dofGID-2;
	k++;

      	values[k] = 0.0;
	colInd[k] = dofGID-1;
	k++;
      }

      if (FDMat_->Filled())
	FDMat_->ReplaceGlobalValues(dofGID, k,
				values.data(), colInd.data());
      else
	FDMat_->InsertGlobalValues(dofGID, k,
			       values.data(), colInd.data());

      // update counter
      l++;
    }//over dof
  }//loop over grid pts

  if(!FDMat_->Filled())
    FDMat_->FillComplete();
}// end method


void UnsteadyNonLinAdvDiffReac2dEpetra::computeChem
( const state_type & C ) const{

  chemForcing_->PutScalar( static_cast<scalar_type>(0) );
  unsigned int k = 0;
  // loop over grid points
  for (auto i=0; i<gptPerProc_; i++){
    // for a given grid point, loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++){
      if (iDof == 0){
  	(*chemForcing_)[k] = -K_*C[k]*C[k+1] + (*s1_)[i];
      }
      else if (iDof == 1){
  	(*chemForcing_)[k] = -K_*C[k-1]*C[k] + (*s2_)[i];
      }
      else if (iDof == 2){
  	(*chemForcing_)[k] = K_*C[k-2]*C[k-1] -K_*C[k] + (*s3_)[i];
      }

      // update indexer
      k++;
    }
  }
}


void UnsteadyNonLinAdvDiffReac2dEpetra::computeJacobian( const state_type & yState ) const
{
  /* the jacobian for this problem can be obtained by
   * computing the FD matrix and then adding the
   * contributions from the chemistry */

  // we need to computes the FD matrix all the time because
  // we use SumIntoGlobalValues below
  this->assembleFDMatrix();

  std::array<int, 3> colInd{};
  std::array<scalar_type, 3> values{};
  int l = {0};
  int numEntries = {0};
  scalar_type c0 = {}, c1 = {};

  for (auto i=0; i<gptPerProc_; i++)
  {
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      numEntries = 0;

      // global ID of current UNKNOWN
      auto dofGID = MyGlobalDof_[l];

      // store the conc values of the species at this point
      if (iDof == 0){
       	c0 = yState[l];	  c1 = yState[l+1];
      }
      if (iDof == 1){
      	c0 = yState[l-1]; c1 = yState[l];
      }
      if (iDof == 2){
      	c0 = yState[l-2]; c1 = yState[l-1];
      }

      // for c_0 we have a -K*c0*c1
      if (iDof == 0){
      	values[0] = -K_*c1; colInd[0] = dofGID;
      	values[1] = -K_*c0; colInd[1] = dofGID+1;
      	numEntries = 2;
      }

      // for c_1, we have -K*c0*c1
      if (iDof == 1){
      	values[0] = -K_*c0; colInd[0] = dofGID;
      	values[1] = -K_*c1; colInd[1] = dofGID-1;
      	numEntries = 2;
      }

      // for c_2, we have K*c0*c1 - K*c3, so three more terms in jacobian
      if (iDof == 2){
      	values[0] = -K_;   colInd[0] = dofGID;
      	values[1] = K_*c1; colInd[1] = dofGID-2;
      	values[2] = K_*c0; colInd[2] = dofGID-1;
      	numEntries = 3;
      }

      FDMat_->SumIntoGlobalValues(dofGID, numEntries,
				  values.data(),
				  colInd.data());

      // update counter
      l++;
    }
  }

  //FDMat_->Print(std::cout);

}//computeJacobian

}} //namespace rompp::apps
