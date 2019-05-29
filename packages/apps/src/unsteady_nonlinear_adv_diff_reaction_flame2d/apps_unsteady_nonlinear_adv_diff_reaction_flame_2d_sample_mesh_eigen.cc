
#include "apps_unsteady_nonlinear_adv_diff_reaction_flame_2d_sample_mesh_eigen.hpp"

namespace rompp{ namespace apps{

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setupPhysicalGrid(){

  // x,y for every mesh cell, this is the set of all cells not just those
  // where we compute residual
  x_.resize(numGpt_);
  y_.resize(numGpt_);

  // u,v are only needed where we compute residual
  u_.resize(numGpt_r_);
  v_.resize(numGpt_r_);
  u_.setConstant(50.0);
  v_.setZero();

  // to tag the region where a cell belongs to
  // this is just needed for left BC
  regionLabel_.resize(numGpt_);

  // the center of the cell at the lower left of the domain,
  // which we use as starting point to create the entire grid
  const auto Ox = this_t::oneHalf*dx_;
  const auto Oy = this_t::oneHalf*dy_;

  // loop over all cells
  for (auto i=0; i<gidsMap_.size(); i++){
    // get the GID of this cell wrt sample mesh
    cellGID_ = gidsMap_[i][0];

    // get the corresponding GID wrt FULL mesh
    cellGIDinFullMesh_ = gidsMap_[i][1];

    // compute the gi, gj of this cell
    globalIDToGiGj(cellGIDinFullMesh_, gi_, gj_);

    // store x,y
    x_[i] = Ox + gi_ * dx_;
    y_[i] = Oy + gj_ * dy_;

    // maybe change to use gj instead of coordinates
    if (y_[i] <= 0.3 ) regionLabel_[i] = -1.0;
    if (y_[i] > .3 and y_[i]<= 0.6) regionLabel_[i] = 0.0;
    if (y_[i] > .6) regionLabel_[i] = 1.0;
  }
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setupFields(){
  //create Jacobian and state
  J_.resize(numDof_r_, numDof_r_);
  state_.resize(numDof_);

  // init to zero
  J_.setZero();
  state_.setZero();

  // the init condition is zero everywhere, except 300 Kelvin for temperature
  // temperature is the first DOF at each cell, so we start from 0 and loop stride of 4
  // to pick the temperature at each cell
  for (auto i=0; i<numDof_; i+=4)
    state_[i] = initTemp_;
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::compute_sources(scalar_type wT,
							    scalar_type wH2,
							    scalar_type wO2,
							    scalar_type wH2O)const{
  const auto H2Ratio	= std::pow(rhoOvWH2 * wH2, 2);
  const auto O2Ratio	= (rhoOvWO2 * wO2);
  const auto expTerm	= preExp_ * std::exp( negE_ /(gasR_ * wT) );
  s_[1]	= -2. * WH2ovRho * H2Ratio * O2Ratio * expTerm;
  s_[2]	= -WO2ovRho * H2Ratio * O2Ratio * expTerm;
  s_[3] = 2. * WH2OovRho * H2Ratio * O2Ratio * expTerm;
  s_[0] = Q_*s_[3];
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::compute_dsdw(scalar_type wT,
								   scalar_type wH2,
								   scalar_type wO2,
								   scalar_type wH2O) const
{
  const auto H2Ratio	     = std::pow(rhoOvWH2 * wH2, this_t::two);
  const auto d_H2Ratio_d_wH2 = this_t::two*(rhoOvWH2 * wH2) * rhoOvWH2;

  const auto O2Ratio	     = (rhoOvWO2 * wO2);
  const auto d_O2Ratio_d_wO2 = rhoOvWO2;

  const auto expTerm	     = preExp_ * std::exp( negE_ /(gasR_ * wT) );
  const auto d_expCoeff_d_wT = -negE_/(gasR_*wT*wT);

  const auto qH2	  = - this_t::two * WH2ovRho * H2Ratio * O2Ratio * expTerm;
  const auto d_qH2_d_wT	  = qH2 * d_expCoeff_d_wT;
  const auto d_qH2_d_wH2  = -this_t::two * WH2ovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
  const auto d_qH2_d_wO2  = -this_t::two * WH2ovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
  const auto d_qH2_d_wH2O = this_t::zero;

  const auto qO2	  = -WO2ovRho * H2Ratio * O2Ratio * expTerm;
  const auto d_qO2_d_wT   = qO2 * d_expCoeff_d_wT;
  const auto d_qO2_d_wH2  = -WO2ovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
  const auto d_qO2_d_wO2  = -WO2ovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
  const auto d_qO2_d_wH2O = this_t::zero;

  const auto qH2O	   = this_t::two * WH2OovRho * H2Ratio * O2Ratio * expTerm;
  const auto d_qH2O_d_wT   = qH2O * d_expCoeff_d_wT;
  const auto d_qH2O_d_wH2  = this_t::two * WH2OovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
  const auto d_qH2O_d_wO2  = this_t::two * WH2OovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
  const auto d_qH2O_d_wH2O = this_t::zero;

  const auto d_qT_d_wT   = Q_ * d_qH2O_d_wT;
  const auto d_qT_d_wH2  = Q_ * d_qH2O_d_wH2;
  const auto d_qT_d_wO2  = Q_ * d_qH2O_d_wO2;
  const auto d_qT_d_wH2O = this_t::zero;

  // first row
  dsdw_(0,0) = d_qT_d_wT;
  dsdw_(0,1) = d_qT_d_wH2;
  dsdw_(0,2) = d_qT_d_wO2;
  dsdw_(0,3) = d_qT_d_wH2O;
  // second row
  dsdw_(1,0) = d_qH2_d_wT;
  dsdw_(1,1) = d_qH2_d_wH2;
  dsdw_(1,2) = d_qH2_d_wO2;
  dsdw_(1,3) = d_qH2_d_wH2O;
  // third row
  dsdw_(2,0) = d_qO2_d_wT;
  dsdw_(2,1) = d_qO2_d_wH2;
  dsdw_(2,2) = d_qO2_d_wO2;
  dsdw_(2,3) = d_qO2_d_wH2O;
  // fourth row
  dsdw_(3,0) = d_qH2O_d_wT;
  dsdw_(3,1) = d_qH2O_d_wH2;
  dsdw_(3,2) = d_qH2O_d_wO2;
  dsdw_(3,3) = d_qH2O_d_wH2O;
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setup(){
  // recall that numGPt_ is the total number of STATE cells,
  // of which those we have residual is a subset
  numGpt_ = gidsMap_.size();
  numDof_ = this_t::numSpecies_ * numGpt_;

  // numGPt_r_ is the number of cells where we want the residual
  // which we can extract from the graph since the graph contains only
  // points where want residual
  numGpt_r_ = graph_.size();
  numDof_r_ = this_t::numSpecies_ * numGpt_r_;

  setupPhysicalGrid();
  setupFields();
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::residual_impl
(const state_type & yState, residual_type & R) const
{
  R.setZero();

  scalar_type c_ip1={}, c_im1={};
  scalar_type c_jp1={}, c_jm1={};
  int dofID = 0;

  // loop over cells where residual needs to be computed
  for (size_t rPt=0; rPt < graph_.size(); ++rPt){

    // global ID of this cell
    const auto cellGID_  = graph_[rPt][0];

    // the gid of this residual grid point seen within the full mesh
    cellGIDinFullMesh_ = gidsMap_[cellGID_][1];

    // find global i,j of this point as if it was in the full mesh
    globalIDToGiGj(cellGIDinFullMesh_, gi_, gj_);

    // the velocity at this cell
    auto uij = u_[cellGID_];
    auto vij = v_[cellGID_];

    // label to help with BC on left boundary
    const auto thisCellLabel = regionLabel_[cellGID_];

    // the state values at current grid point
    const auto wT   = yState[cellGID_*numSpecies_];
    const auto wH2  = yState[cellGID_*numSpecies_+1];
    const auto wO2  = yState[cellGID_*numSpecies_+2];
    const auto wH2O = yState[cellGID_*numSpecies_+3];
    // compute source terms at this location
    compute_sources(wT, wH2, wO2, wH2O);


    // get ID of the first state dof (temperature) for the cell on the left
    const auto firstDofID_westCell  = graph_[rPt][1]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell above
    const auto firstDofID_northCell = graph_[rPt][2]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell on right
    const auto firstDofID_eastCell  = graph_[rPt][3]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell below
    const auto firstDofID_southCell = graph_[rPt][4]*numSpecies_;
    // so to get the other dofs at each neighboring cell, we just increment
    // for the second state dof (fraction of H2 ) just add 1 to any of the above
    // for the third state dof (fraction of O2 ) just add 2 to any of the above
    // for the fourth state dof (fraction of H2O ) just add 3 to any of the above


    // loop over local dofs and compute residual
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      // this is the global ID of the currrent DOF at the current cell
      dofID = cellGID_ * numSpecies_ + iDof;

      // add diagonal
      R[dofID] = -this_t::two*(K_ * dxSqInv_ + K_ * dySqInv_) * yState[dofID];

      //---------------
      // x-direction
      //---------------
      // near left wall
      if (gi_==0){

	c_ip1 = yState[firstDofID_eastCell + iDof];

	// contribution due to stencil at BC
      	const auto leftBC = (thisCellLabel != 0) ?
      	  bcLeftGamma13_[iDof] : bcLeftGamma2_[iDof];
      	R[dofID] += -this_t::two*(K_*dxSqInv_ + uij*dx2Inv_)*yState[dofID];
      	R[dofID] += this_t::oneThird*(K_*dxSqInv_ + uij*dx2Inv_)*c_ip1;
      	R[dofID] += this_t::eightOvThree*(K_*dxSqInv_ + uij*dx2Inv_)*leftBC;

	// contribution from adjecent cell on the right at i+1,j
      	R[dofID] += (K_*dxSqInv_ - uij*dx2Inv_)*c_ip1;
      }

      // i-1,j and i+1,j (interior points, regular stencil)
      if (gi_>=1 and gi_<Nx_-1){

	// contribution from left cell
      	c_im1 = yState[firstDofID_westCell + iDof];
      	R[dofID] += (K_*dxSqInv_ + uij*dx2Inv_)*c_im1;

	// contribution from right cell
	c_ip1 = yState[firstDofID_eastCell + iDof];
      	R[dofID] += (K_*dxSqInv_ - uij*dx2Inv_)*c_ip1;
      }

      // right wall we have homog Neumann BC
      if (gi_==Nx_-1){
	c_im1 = yState[firstDofID_westCell + iDof];
      	R[dofID] += this_t::two*K_*dxSqInv_ * c_im1;
      }

      //---------------
      // y-direction
      //---------------
      // i, j-1 and i, j+1
      if (gj_>=1 and gj_<Ny_-1){
	// contribution from cell below
      	c_jm1 = yState[firstDofID_southCell + iDof];
      	R[dofID] += (K_*dySqInv_ + vij*dy2Inv_)*c_jm1;

	// contribution from cell above
      	c_jp1 = yState[firstDofID_northCell + iDof];
      	R[dofID] += (K_*dySqInv_ - vij*dy2Inv_)*c_jp1;
      }

      // bottom wall we have homog Neumann BC
      if (gj_==0){
      	c_jp1 = yState[firstDofID_northCell + iDof];
      	R[dofID] += this_t::two*K_*dySqInv_ * c_jp1;
      }

      // top wall we have homog Neumann BC
      if (gj_==Ny_-1){
      	c_jm1 = yState[firstDofID_southCell + iDof];
      	R[dofID] += this_t::two*K_*dySqInv_ * c_jm1;
      }

      // account for sources
      R[dofID] += s_[iDof];

    }//loop over dof
  }//loop over grid pts
}// end method



void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::jacobian_impl
(const state_type & yState, jacobian_type & jac)const{

  if (jac.rows() == 0 || jac.cols()==0 )
    jac.resize(numDof_r_, numDof_r_);

  // triplets is used to store a series of (row, col, value)
  // has to be cleared becuase we append to it while computing
  tripletList.clear();

  scalar_type value = {0};
  int dofID = 0;

  // loop over cells where residual needs to be computed
  for (size_t rPt=0; rPt < graph_.size(); ++rPt){

    // global ID of this cell
    const auto cellGID_  = graph_[rPt][0];

    // the gid of this residual grid point seen within the full mesh
    cellGIDinFullMesh_ = gidsMap_[cellGID_][1];

    // find global i,j of this point as if it was in the full mesh
    globalIDToGiGj(cellGIDinFullMesh_, gi_, gj_);

    // the velocity at this cell
    auto uij = u_[cellGID_];
    auto vij = v_[cellGID_];

    // the state at current grid point
    // the state values at current grid point
    const auto wT   = yState[cellGID_*numSpecies_];
    const auto wH2  = yState[cellGID_*numSpecies_+1];
    const auto wO2  = yState[cellGID_*numSpecies_+2];
    const auto wH2O = yState[cellGID_*numSpecies_+3];
    compute_dsdw(wT, wH2, wO2, wH2O);

    // get ID of the first state dof (temperature) for the cell on the left
    const auto firstDofID_westCell  = graph_[rPt][1]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell above
    const auto firstDofID_northCell = graph_[rPt][2]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell on right
    const auto firstDofID_eastCell  = graph_[rPt][3]*numSpecies_;
    // get ID of the first state dof (temperature) for the cell below
    const auto firstDofID_southCell = graph_[rPt][4]*numSpecies_;
    // so to get the other dofs at each neighboring cell, we just increment
    // for the second state dof (fraction of H2 ) just add 1 to any of the above
    // for the third state dof (fraction of O2 ) just add 2 to any of the above
    // for the fourth state dof (fraction of H2O ) just add 3 to any of the above


    // loop over local dofs
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      // this is the global ID of the currrent DOF at the current cell
      dofID = cellGID_ * numSpecies_ + iDof;

      // diagonal
      auto diagVal = -this_t::two*(K_ * dxSqInv_ + K_ * dySqInv_) + dsdw_(iDof,iDof);

      // if we are near the left wall we need to account from contributions
      // stemming from how we treat the BC
      if (gi_==0){
	// contributions coming from how to treat left BC
	diagVal += -this_t::two*(K_*dxSqInv_ + uij*dx2Inv_);

      	value = this_t::oneThird*(K_*dxSqInv_ + uij*dx2Inv_);
	// jacobian contribution from left BC treatment relative to cell on the right
      	tripletList.push_back( Tr( dofID, firstDofID_eastCell+iDof, value) );
      }
      tripletList.push_back( Tr( dofID, dofID, diagVal) );

      // i-1, j and i+1, j
      if (gi_>=1 and gi_<Nx_-1){
	// jacob wrt left cell
	value = (K_*dxSqInv_ + uij*dx2Inv_);
	tripletList.push_back( Tr( dofID, firstDofID_westCell + iDof, value) );

	// jacob wrt right cell
	value = (K_*dxSqInv_ - uij*dx2Inv_);
	tripletList.push_back( Tr( dofID, firstDofID_eastCell + iDof, value) );
      }

      // right wall we have homog Neumann BC
      if (gi_==Nx_-1){
      	value = this_t::two*K_*dxSqInv_;
	tripletList.push_back( Tr( dofID, firstDofID_westCell + iDof, value) );
      }

      // i, j-1 and i, j+1
      if (gj_>=1 and gj_<Ny_-1){
	value = (K_*dySqInv_ + vij*dy2Inv_);
	tripletList.push_back( Tr( dofID, firstDofID_southCell + iDof, value) );

	value = (K_*dySqInv_ - vij*dy2Inv_);
	tripletList.push_back( Tr( dofID, firstDofID_northCell + iDof, value) );
      }

      // bottom wall we have homog Neumann BC
      if (gj_==0){
      	value = this_t::two*K_*dySqInv_;
      	tripletList.push_back( Tr( dofID, firstDofID_northCell + iDof, value) );
      }

      // top wall we have homog Neumann BC
      if (gj_==Ny_-1){
      	value = this_t::two*K_*dySqInv_;
      	tripletList.push_back( Tr( dofID, firstDofID_southCell + iDof, value) );
      }

      // sources contribution
      if (iDof==0){
	tripletList.push_back( Tr( dofID, dofID+1, dsdw_(iDof,1)) );
	tripletList.push_back( Tr( dofID, dofID+2, dsdw_(iDof,2)) );
	tripletList.push_back( Tr( dofID, dofID+3, dsdw_(iDof,3)) );
      }
      if (iDof==1){
	tripletList.push_back( Tr( dofID, dofID-1, dsdw_(iDof,0)) );
	tripletList.push_back( Tr( dofID, dofID+1, dsdw_(iDof,2)) );
	tripletList.push_back( Tr( dofID, dofID+2, dsdw_(iDof,3)) );
      }
      if (iDof == 2){
	tripletList.push_back( Tr( dofID, dofID-2, dsdw_(iDof,0)) );
	tripletList.push_back( Tr( dofID, dofID-1, dsdw_(iDof,1)) );
	tripletList.push_back( Tr( dofID, dofID+1, dsdw_(iDof,3)) );
      }
      if (iDof == 3){
	tripletList.push_back( Tr( dofID, dofID-3, dsdw_(iDof,0)) );
	tripletList.push_back( Tr( dofID, dofID-2, dsdw_(iDof,1)) );
	tripletList.push_back( Tr( dofID, dofID-1, dsdw_(iDof,2)) );
      }

    }// dof loop
  }

  // use triplets to fill the jacobian
  jac.setFromTriplets(tripletList.begin(), tripletList.end());

}//end jacob

}} //namespace rompp::apps
