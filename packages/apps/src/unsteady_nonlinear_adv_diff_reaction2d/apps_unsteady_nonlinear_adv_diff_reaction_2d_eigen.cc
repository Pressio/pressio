
#include "apps_unsteady_nonlinear_adv_diff_reaction_2d_eigen.hpp"

namespace rompp{ namespace apps{

void UnsteadyNonLinAdvDiffReac2dEigen::setupPhysicalGrid(){
  x_.resize(numGpt_);
  y_.resize(numGpt_);
  u_.resize(numGpt_);
  v_.resize(numGpt_);

  int gi{}; int gj{};
  for (auto i = 0; i<numGpt_; i++){
    auto GID = i;
    globalIDToGiGj(GID, gi, gj);
    x_[i] = dx_ + gi * dx_;
    y_[i] = gj * dy_;
    auto xval = x_[i];
    auto yval = y_[i];
    u_[i] = -std::sin(M_PI*xval) * std::cos(M_PI*yval);
    v_[i] = std::cos(M_PI*xval) * std::sin(M_PI*yval);
  }
}

void UnsteadyNonLinAdvDiffReac2dEigen::setupFields(){
  s1_.resize(numGpt_);
  s2_.resize(numGpt_);
  s3_.resize(numGpt_);
  J_.resize(numDof_, numDof_);
  state_.resize(numDof_);
  J_.setZero();
  state_.setZero();
}

void UnsteadyNonLinAdvDiffReac2dEigen::fillSource1(){
  for (auto i=0; i<numGpt_; i++){
    const scalar_type xij = x_[i];
    const scalar_type yij = y_[i];
    const scalar_type dx = (xij - oPtS1[0]);
    const scalar_type dy = (yij - oPtS1[1]);
    const scalar_type distance = dx*dx + dy*dy;
    s1_[i] = ( std::sqrt(distance) <= rS1) ? 0.1 : 0.0;
  }
}

void UnsteadyNonLinAdvDiffReac2dEigen::fillSource2(){
  for (auto i=0; i<numGpt_; i++){
    const scalar_type xij = x_[i];
    const scalar_type yij = y_[i];
    const scalar_type dx = (xij - oPtS2[0]);
    const scalar_type dy = (yij - oPtS2[1]);
    const scalar_type distance = dx*dx + dy*dy;
    s2_[i] = ( std::sqrt(distance) <= rS2) ? 0.1 : 0.0;
  }
}

void UnsteadyNonLinAdvDiffReac2dEigen::fillSource3(){
  s3_.setZero();
}

void UnsteadyNonLinAdvDiffReac2dEigen::setup(){
  numGpt_ = Nx_ * Ny_;
  numDof_ = this_t::numSpecies_ * Nx_ * Ny_;
  setupPhysicalGrid();
  setupFields();
  fillSource1();
  fillSource2();
}

void UnsteadyNonLinAdvDiffReac2dEigen::residual_impl
(const state_type & yState, residual_type & R) const
{
  /* note that R is the residual vector, which has
   * a dofMap_ because it contains all dofs.
   * Once again, the dofs are the numFields * numOfUnkownGridPoints
   */
  int gi{}; int gj{};

  scalar_type c_ip1={}, c_im1={};
  scalar_type c_jp1={}, c_jm1={};

  // global ID of current UNKNOWN
  int dofGID = {0};

  // rename diffusivity for convenience
  auto D = eps_;

  // loop over grid points
  for (auto GID=0; GID<numGpt_; ++GID){
    // the local valocity
    auto uij = u_[GID];
    auto vij = v_[GID];

    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    // for a given grid point, loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      // add diagonal
      R[dofGID] = -2.0*(D * dxSqInv_ + D * dySqInv_) * yState[dofGID];

      // i-1, j
      if (gi>=1){
      	c_im1 = yState[dofGID-numSpecies_];
      	R[dofGID] += (D*dxSqInv_ + uij*dx2Inv_)*c_im1;
      }

      // i+1, j
      if (gi<Nx_-1){
      	c_ip1 = yState[dofGID+numSpecies_];
      	R[dofGID] += (D*dxSqInv_ - uij*dx2Inv_)*c_ip1;
      }

      // i, j-1 and i, j+1
      if (gj>=1 and gj<Ny_-1){
      	c_jm1 = yState[dofGID - Nx_*numSpecies_];
      	R[dofGID] += (D*dySqInv_ + vij*dy2Inv_)*c_jm1;

      	c_jp1 = yState[dofGID + Nx_*numSpecies_];
      	R[dofGID] += (D*dySqInv_ - vij*dy2Inv_)*c_jp1;
      }

      // bottom wall we have homog Neumann BC
      if (gj==0){
      	c_jp1 = yState[dofGID + Nx_*numSpecies_];
      	R[dofGID] += 2.0*D*dySqInv_ * c_jp1;
      }

      // top wall we have homog Neumann BC
      if (gj==Ny_-1){
      	c_jm1 = yState[dofGID - Nx_*numSpecies_];
      	R[dofGID] += 2.0*D*dySqInv_ * c_jm1;
      }

      // account for chemistry and sources
      if (iDof == 0)
      	R[dofGID] += -K_*yState[dofGID]*yState[dofGID+1] + s1_[GID];

      if (iDof == 1)
      	R[dofGID] += -K_*yState[dofGID-1]*yState[dofGID] + s2_[GID];

      if (iDof == 2){
      	R[dofGID] += K_*yState[dofGID-2]*yState[dofGID-1]
      	  - K_* yState[dofGID] + s3_[GID];
      }

      // update counter
      dofGID++;

    }//loop over dof
  }//loop over grid pts

}// end method



void UnsteadyNonLinAdvDiffReac2dEigen::jacobian_impl
(const state_type & yState, jacobian_type & jac)const{

  if (jac.rows() == 0 || jac.cols()==0 ){
    jac.resize(yState.size(), yState.size());
  }
  // triplets is just to store a series of (row,col,value)
  tripletList.clear();

  int gi{}; int gj{};
  int dofGID = {0};
  auto D = eps_;

  scalar_type value = {0};

  // loop over grid points
  for (auto GID=0; GID<numGpt_; ++GID){
    auto uij = u_[GID]; auto vij = v_[GID];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    // loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      // add diagonal
      auto value1 = -2.0*(D * dxSqInv_ + D * dySqInv_);
      scalar_type value2{0};
      if (iDof==0) value2 = -K_*yState[dofGID+1];
      if (iDof==1) value2 = -K_*yState[dofGID-1];
      if (iDof==2) value2 = -K_;
      tripletList.push_back( Tr( dofGID, dofGID, value1+value2) );

      // i-1, j
      if (gi>=1){
	value = (D*dxSqInv_ + uij*dx2Inv_);
	auto im1_col = dofGID-numSpecies_;
	tripletList.push_back( Tr( dofGID, im1_col, value) );
      }

      // i+1, j
      if (gi<Nx_-1){
	value = (D*dxSqInv_ - uij*dx2Inv_);
	auto ip1_col = dofGID+numSpecies_;
	tripletList.push_back( Tr( dofGID, ip1_col, value) );
      }

      // i, j-1 and i, j+1
      if (gj>=1 and gj<Ny_-1){
	value = (D*dySqInv_ + vij*dy2Inv_);
	auto jm1_col = dofGID - Nx_*numSpecies_;
	tripletList.push_back( Tr( dofGID, jm1_col, value) );

	value = (D*dySqInv_ - vij*dy2Inv_);
	auto jp1_col = dofGID + Nx_*numSpecies_;
	tripletList.push_back( Tr( dofGID, jp1_col, value) );
      }

      // bottom wall we have homog Neumann BC
      if (gj==0){
      	value = 2.0*D*dySqInv_;
      	auto jp1_col = dofGID + Nx_*numSpecies_;
      	tripletList.push_back( Tr( dofGID, jp1_col, value) );
      }

      // top wall we have homog Neumann BC
      if (gj==Ny_-1){
      	value = 2.0*D*dySqInv_;
      	auto jm1_col = dofGID - Nx_*numSpecies_;
      	tripletList.push_back( Tr( dofGID, jm1_col, value) );
      }

      // account for chemistry (products of conc for dof = 2)
      if (iDof==0){
	value = -K_*yState[dofGID];
	tripletList.push_back( Tr( dofGID, dofGID+1, value) );
      }
      if (iDof==1){
	value = -K_*yState[dofGID];
	tripletList.push_back( Tr( dofGID, dofGID-1, value) );
      }
      if (iDof == 2){
      	value = K_ * yState[dofGID-1];
      	tripletList.push_back( Tr( dofGID, dofGID-2, value) );

      	value = K_ * yState[dofGID-2];
      	tripletList.push_back( Tr( dofGID, dofGID-1, value) );
      }

      dofGID++;
    }// dof loop
  }

  // use triplets to fill the jacobian
  jac.setFromTriplets(tripletList.begin(), tripletList.end());

}//end jacob

}} //namespace rompp::apps
