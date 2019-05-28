
#include "apps_unsteady_nonlinear_adv_diff_reaction_flame_2d_sample_mesh_eigen.hpp"

namespace rompp{ namespace apps{

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setupPhysicalGrid(){
  x_.resize(numGpt_);
  y_.resize(numGpt_);
  u_.resize(numGpt_);
  v_.resize(numGpt_);
  regionLabel_.resize(numGpt_);

  // the grid point at the lower left
  const auto Ox = 0.5*dx_;
  const auto Oy = 0.5*dy_;

  int gi{}; int gj{};
  int GID = 0;
  for (auto i=0; i<gidsMap_.size(); i++){
    GID = gidsMap_[i][0];
    auto GIDinFullMesh = gidsMap_[i][1];
    globalIDToGiGj(GIDinFullMesh, gi, gj);

    x_[i] = Ox + gi * dx_;
    y_[i] = Oy + gj * dy_;
    u_[i] = 50.0;
    v_[i] = 0.0;

    // maybe change to use gj instead of coordinates
    if (y_[i] <= 0.3 ) regionLabel_[i] = -1.0;
    if (y_[i] > .3 and y_[i]<= 0.6) regionLabel_[i] = 0.0;
    if (y_[i] > .6) regionLabel_[i] = 1.0;
  }
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setupFields(){
  //J_.resize(numDof_, numDof_);
  //J_.setZero();
  state_.resize(numDof_);

  // the init condition is zero everywhere, but 300 Kelvin for temperature
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

// void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::compute_dsdw(scalar_type wT,
// 							 scalar_type wH2,
// 							 scalar_type wO2,
// 							 scalar_type wH2O)const{
//   const auto H2Ratio	     = std::pow(rhoOvWH2 * wH2, 2);
//   const auto d_H2Ratio_d_wH2 = 2*(rhoOvWH2 * wH2) * rhoOvWH2;

//   const auto O2Ratio	     = (rhoOvWO2 * wO2);
//   const auto d_O2Ratio_d_wO2 = rhoOvWO2;

//   const auto expTerm	     = preExp_ * std::exp( negE_ /(gasR_ * wT) );
//   const auto d_expCoeff_d_wT = -negE_/(gasR_*wT*wT);

//   const auto qH2	  = -2. * WH2ovRho * H2Ratio * O2Ratio * expTerm;
//   const auto d_qH2_d_wT	  = qH2 * d_expCoeff_d_wT;
//   const auto d_qH2_d_wH2  = -2. * WH2ovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
//   const auto d_qH2_d_wO2  = -2. * WH2ovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
//   const auto d_qH2_d_wH2O = 0.0;

//   const auto qO2	  = -WO2ovRho * H2Ratio * O2Ratio * expTerm;
//   const auto d_qO2_d_wT   = qO2 * d_expCoeff_d_wT;
//   const auto d_qO2_d_wH2  = -WO2ovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
//   const auto d_qO2_d_wO2  = -WO2ovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
//   const auto d_qO2_d_wH2O = 0.0;

//   const auto qH2O	   = 2. * WH2OovRho * H2Ratio * O2Ratio * expTerm;
//   const auto d_qH2O_d_wT   = qH2O * d_expCoeff_d_wT;
//   const auto d_qH2O_d_wH2  = 2. * WH2OovRho * d_H2Ratio_d_wH2 * O2Ratio * expTerm;
//   const auto d_qH2O_d_wO2  = 2. * WH2OovRho * H2Ratio * d_O2Ratio_d_wO2 * expTerm;
//   const auto d_qH2O_d_wH2O = 0.0;

//   const auto d_qT_d_wT   = Q_ * d_qH2O_d_wT;
//   const auto d_qT_d_wH2  = Q_ * d_qH2O_d_wH2;
//   const auto d_qT_d_wO2  = Q_ * d_qH2O_d_wO2;
//   const auto d_qT_d_wH2O = 0.0;

//   dsdw_(0,0) = d_qT_d_wT;
//   dsdw_(0,1) = d_qT_d_wH2;
//   dsdw_(0,2) = d_qT_d_wO2;
//   dsdw_(0,3) = d_qT_d_wH2O;

//   dsdw_(1,0) = d_qH2_d_wT;
//   dsdw_(1,1) = d_qH2_d_wH2;
//   dsdw_(1,2) = d_qH2_d_wO2;
//   dsdw_(1,3) = d_qH2_d_wH2O;

//   dsdw_(2,0) = d_qO2_d_wT;
//   dsdw_(2,1) = d_qO2_d_wH2;
//   dsdw_(2,2) = d_qO2_d_wO2;
//   dsdw_(2,3) = d_qO2_d_wH2O;

//   dsdw_(3,0) = d_qH2O_d_wT;
//   dsdw_(3,1) = d_qH2O_d_wH2;
//   dsdw_(3,2) = d_qH2O_d_wO2;
//   dsdw_(3,3) = d_qH2O_d_wH2O;
// }

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::setup(){
  numGpt_ = gidsMap_.size();
  numDof_ = this_t::numSpecies_ * numGpt_;

  numGpt_r_ = graph_.size();
  numDof_r_ = this_t::numSpecies_ * numGpt_r_;

  setupPhysicalGrid();
  setupFields();
}

void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::residual_impl
(const state_type & yState, residual_type & R) const
{
  /* note that R is the residual vector, which has
   * a dofMap_ because it contains all dofs.
   * Once again, the dofs are the numFields * numOfUnkownGridPoints
   */
  R.setZero();

  int gi{}; int gj{};
  scalar_type c_ip1={}, c_im1={};
  scalar_type c_jp1={}, c_jm1={};

  // diffusivity alias for convenience
  const auto D = K_;

  // loop over grid points where residual needs to be computed
  for (size_t rPt=0; rPt < graph_.size(); ++rPt){

    // get from the graph the gID of this point
    auto rPtGID  = graph_[rPt][0];

    int dofGID = rPtGID*numSpecies_;

    // the id of this residual grid point seen within the full mesh
    auto rPt_fmesh = gidsMap_[rPtGID][1];

    // the velocity at this point
    auto uij = u_[rPtGID];
    auto vij = v_[rPtGID];

    // label to help with BC on left boundary
    const auto myLabel = regionLabel_[rPtGID];

    // find global i,j of this point as if it was in the full mesh
    globalIDToGiGj(rPt_fmesh, gi, gj);

    // the state values at current grid point
    const auto wT   = yState[dofGID];
    const auto wH2  = yState[dofGID+1];
    const auto wO2  = yState[dofGID+2];
    const auto wH2O = yState[dofGID+3];
    // compute source terms at this location
    compute_sources(wT, wH2, wO2, wH2O);

    // get ID of left point
    auto dofGID_east  = graph_[rPt][1]*numSpecies_;
    // get ID of top point
    auto dofGID_north  = graph_[rPt][2]*numSpecies_;
    // get ID of right point
    auto dofGID_west  = graph_[rPt][3]*numSpecies_;
    // get ID of bottom point
    auto dofGID_south  = graph_[rPt][4]*numSpecies_;

    // loop over local dofs and compute residual
    for (auto iDof=0; iDof<numSpecies_; iDof++)
    {
      // add diagonal
      R[dofGID] = -2.0*(D * dxSqInv_ + D * dySqInv_) * yState[dofGID];

      // near left wall
      if (gi==0){

      	const auto leftBC = (myLabel != 0) ?
      	  bcLeftGamma13_[iDof] : bcLeftGamma2_[iDof];

      	R[dofGID] += -2.0*(D*dxSqInv_ + uij*dx2Inv_)*yState[dofGID];
      	R[dofGID] += (1./3.)*(D*dxSqInv_ + uij*dx2Inv_)*yState[dofGID+numSpecies_];
      	R[dofGID] += (8./3.)*(D*dxSqInv_ + uij*dx2Inv_)*leftBC;

      	c_ip1 = yState[dofGID + numSpecies_];
      	R[dofGID] += (D*dxSqInv_ - uij*dx2Inv_)*c_ip1;
      }

      // i-1,j and i+1,j
      if (gi>=1 and gi<Nx_-1){
      	c_im1 = yState[dofGID-numSpecies_];
      	R[dofGID] += (D*dxSqInv_ + uij*dx2Inv_)*c_im1;

      	c_ip1 = yState[dofGID+numSpecies_];
      	R[dofGID] += (D*dxSqInv_ - uij*dx2Inv_)*c_ip1;
      }

      // right wall we have homog Neumann BC
      if (gi==Nx_-1){
      	c_im1 = yState[dofGID - numSpecies_];
      	R[dofGID] += 2.0*D*dxSqInv_ * c_im1;
      }

      // i, j-1 and i, j+1
      if (gj>=1 and gj<Ny_-1){
      	c_jm1 = yState[dofGID_south + iDof];
      	R[dofGID] += (D*dySqInv_ + vij*dy2Inv_)*c_jm1;

      	c_jp1 = yState[dofGID_north + iDof];
      	R[dofGID] += (D*dySqInv_ - vij*dy2Inv_)*c_jp1;
      }

      // bottom wall we have homog Neumann BC
      if (gj==0){
      	c_jp1 = yState[dofGID_north + iDof];
      	R[dofGID] += 2.0*D*dySqInv_ * c_jp1;
      }

      // top wall we have homog Neumann BC
      if (gj==Ny_-1){
      	c_jm1 = yState[dofGID_south + iDof];
      	R[dofGID] += 2.0*D*dySqInv_ * c_jm1;
      }

      // account for sources
      R[dofGID] += s_[iDof];

      // update counter
      dofGID++;

    }//loop over dof
  }//loop over grid pts
}// end method



// void UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen::jacobian_impl
// (const state_type & yState, jacobian_type & jac)const{

//   if (jac.rows() == 0 || jac.cols()==0 ){
//     jac.resize(yState.size(), yState.size());
//   }
//   // triplets is just to store a series of (row,col,value)
//   tripletList.clear();

//   int gi{}; int gj{};
//   int dofGID = {0};
//   auto D = K_;

//   scalar_type value = {0};

//   // loop over grid points
//   for (auto GID=0; GID<numGpt_; ++GID){
//     auto uij = u_[GID];
//     auto vij = v_[GID];

//     // find global i,j of this point
//     globalIDToGiGj(GID, gi, gj);

//     // the state at current grid point
//     const auto wT   = yState[dofGID];
//     const auto wH2  = yState[dofGID+1];
//     const auto wO2  = yState[dofGID+2];
//     const auto wH2O = yState[dofGID+3];
//     compute_dsdw(wT, wH2, wO2, wH2O);

//     // loop over local dofs
//     for (auto iDof=0; iDof<numSpecies_; iDof++)
//     {
//       // diagonal
//       auto value1 = -2.0*(D * dxSqInv_ + D * dySqInv_);
//       auto diagVal = value1 + dsdw_(iDof,iDof);

//       if (gi==0){
// 	diagVal -= 2.0*(D*dxSqInv_ + uij*dx2Inv_);

//       	value = (D*dxSqInv_ + uij*dx2Inv_)*(1./3.);
//       	auto ip1_col = dofGID+numSpecies_;
//       	tripletList.push_back( Tr( dofGID, ip1_col, value) );
//       }
//       tripletList.push_back( Tr( dofGID, dofGID, diagVal) );


//       // i-1, j and i+1, j
//       if (gi>=1 and gi<Nx_-1){
// 	value = (D*dxSqInv_ + uij*dx2Inv_);
// 	auto im1_col = dofGID-numSpecies_;
// 	tripletList.push_back( Tr( dofGID, im1_col, value) );

// 	value = (D*dxSqInv_ - uij*dx2Inv_);
// 	auto ip1_col = dofGID+numSpecies_;
// 	tripletList.push_back( Tr( dofGID, ip1_col, value) );
//       }

//       // right wall
//       if (gi==Nx_-1){
//       	value = 2*D*dxSqInv_;
// 	auto im1_col = dofGID-numSpecies_;
// 	tripletList.push_back( Tr( dofGID, im1_col, value) );
//       }

//       // i, j-1 and i, j+1
//       if (gj>=1 and gj<Ny_-1){
// 	value = (D*dySqInv_ + vij*dy2Inv_);
// 	auto jm1_col = dofGID - Nx_*numSpecies_;
// 	tripletList.push_back( Tr( dofGID, jm1_col, value) );

// 	value = (D*dySqInv_ - vij*dy2Inv_);
// 	auto jp1_col = dofGID + Nx_*numSpecies_;
// 	tripletList.push_back( Tr( dofGID, jp1_col, value) );
//       }

//       // bottom wall we have homog Neumann BC
//       if (gj==0){
//       	value = 2.0*D*dySqInv_;
//       	auto jp1_col = dofGID + Nx_*numSpecies_;
//       	tripletList.push_back( Tr( dofGID, jp1_col, value) );
//       }

//       // top wall we have homog Neumann BC
//       if (gj==Ny_-1){
//       	value = 2.0*D*dySqInv_;
//       	auto jm1_col = dofGID - Nx_*numSpecies_;
//       	tripletList.push_back( Tr( dofGID, jm1_col, value) );
//       }

//       // sources contribution
//       if (iDof==0){
// 	tripletList.push_back( Tr( dofGID, dofGID+1, dsdw_(iDof,1)) );
// 	tripletList.push_back( Tr( dofGID, dofGID+2, dsdw_(iDof,2)) );
// 	tripletList.push_back( Tr( dofGID, dofGID+3, dsdw_(iDof,3)) );
//       }
//       if (iDof==1){
// 	tripletList.push_back( Tr( dofGID, dofGID-1, dsdw_(iDof,0)) );
// 	tripletList.push_back( Tr( dofGID, dofGID+1, dsdw_(iDof,2)) );
// 	tripletList.push_back( Tr( dofGID, dofGID+2, dsdw_(iDof,3)) );
//       }
//       if (iDof == 2){
// 	tripletList.push_back( Tr( dofGID, dofGID-2, dsdw_(iDof,0)) );
// 	tripletList.push_back( Tr( dofGID, dofGID-1, dsdw_(iDof,1)) );
// 	tripletList.push_back( Tr( dofGID, dofGID+1, dsdw_(iDof,3)) );
//       }
//       if (iDof == 3){
// 	tripletList.push_back( Tr( dofGID, dofGID-3, dsdw_(iDof,0)) );
// 	tripletList.push_back( Tr( dofGID, dofGID-2, dsdw_(iDof,1)) );
// 	tripletList.push_back( Tr( dofGID, dofGID-1, dsdw_(iDof,2)) );
//       }

//       dofGID++;
//     }// dof loop
//   }

//   // use triplets to fill the jacobian
//   jac.setFromTriplets(tripletList.begin(), tripletList.end());

// }//end jacob

}} //namespace rompp::apps
