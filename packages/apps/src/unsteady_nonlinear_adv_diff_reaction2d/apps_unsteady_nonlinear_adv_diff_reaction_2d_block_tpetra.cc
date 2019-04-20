
#include "apps_unsteady_nonlinear_adv_diff_reaction_2d_block_tpetra.hpp"

namespace rompp{ namespace apps{

void UnsteadyNonLinAdvDiffReac2dBlockTpetra::createMap(){
  // total number of unknown grid points (we only consider the interior points)
  numGlobalGpt_ = Nx_ * Ny_;
  map_ = Teuchos::rcp(new map_t(numGlobalGpt_, 0, comm_));

  // how many grid points I own
  gptPerProc_ = map_->getNodeNumElements();
  const auto myMinGId = map_->getMinGlobalIndex();
  myGIDs_.resize(gptPerProc_);
  std::iota(myGIDs_.begin(), myGIDs_.end(), myMinGId);
}


void UnsteadyNonLinAdvDiffReac2dBlockTpetra::initGridAndVel(){
  x_ = std::make_shared<tpetVec>(map_);
  y_ = std::make_shared<tpetVec>(map_);
  u_ = std::make_shared<tpetVec>(map_);
  v_ = std::make_shared<tpetVec>(map_);

  auto x_ta = x_->getDataNonConst();
  auto y_ta = y_->getDataNonConst();
  auto u_ta = u_->getDataNonConst();
  auto v_ta = v_->getDataNonConst();
  int gi{}; int gj{};
  for (auto i = 0; i<gptPerProc_; i++){
    auto GID = myGIDs_[i];
    globalIDToGiGj(GID, gi, gj);
    x_ta[i] = dx_ + gi * dx_;
    y_ta[i] = gj * dy_;
    auto xval = x_ta[i];
    auto yval = y_ta[i];
    u_ta[i] = -std::sin(M_PI*xval) * std::cos(M_PI*yval);
    v_ta[i] = std::cos(M_PI*xval) * std::sin(M_PI*yval);
  }
}


void UnsteadyNonLinAdvDiffReac2dBlockTpetra::createGraph(){
  constexpr auto nnZ = UnsteadyNonLinAdvDiffReac2dBlockTpetra::maxNonZeroPerRow_;

  // this is the graph of the grid
  graph_ = std::make_shared<graph_t>(map_, nnZ);

  // my integer "coordinates" in the global 2d grid
  int gi{}; int gj{};

  Teuchos::Array<GO> colInd;
  for (auto i=0; i<gptPerProc_; i++)
  {
    // global ID of current point
    auto GID = myGIDs_[i];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    colInd.push_back(GID);
    if (gi>=1)
      colInd.push_back(GID - 1);

    if (gi<Nx_-1)
      colInd.push_back(GID+1);

    if (gj>=1 and gj<Ny_-1){
      colInd.push_back(GID - Nx_);
      colInd.push_back(GID + Nx_);
    }

    if (gj==0)
      colInd.push_back(GID + Nx_);

    if (gj==Ny_-1)
      colInd.push_back(GID - Nx_);

    graph_->insertGlobalIndices(GID, colInd());
    colInd.clear();

    /* std::cout << rank_ << " row=" << GID << " " << k << " "; */
    /* for (auto l = 0; l < k; ++l) */
    /*   std::cout << colInd[l] << " "; */
    /* std::cout << std::endl; */
  }// rows loop

  graph_->fillComplete();
  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  graph_->describe(*out, Teuchos::VERB_EXTREME);
}


void UnsteadyNonLinAdvDiffReac2dBlockTpetra::initFields(){
  constexpr auto nDof = UnsteadyNonLinAdvDiffReac2dBlockTpetra::numSpecies_;
  constexpr auto zero = static_cast<ST>(0);

  A_ = std::make_shared<nativeMatrix>( *graph_, nDof);
  chemReac_ = std::make_shared<nativeVec>(*map_, nDof);
  state_ = std::make_shared<nativeVec>(*map_, nDof);
  source_ = std::make_shared<nativeVec>(*map_, nDof);

  // set all the zero
  chemReac_->putScalar(zero);
  state_->putScalar(zero);
  source_->putScalar(zero);
}


void UnsteadyNonLinAdvDiffReac2dBlockTpetra::computeSource(){

  std::array<ST, UnsteadyNonLinAdvDiffReac2dBlockTpetra::numSpecies_> values{};
  const auto xta = x_->getData();
  const auto yta = y_->getData();

  scalar_type xij  = {};
  scalar_type yij  = {};
  scalar_type dx   = {};
  scalar_type dy   = {};
  scalar_type dist = {};

  for (auto i=0; i<gptPerProc_; i++){
    xij = xta[i];
    yij = yta[i];

    dx = (xij - oPtS1[0]);
    dy = (yij - oPtS1[1]);
    dist = dx*dx + dy*dy;
    values[0] = ( std::sqrt(dist) <= rS1) ? 0.1 : 0.0;

    dx = (xij - oPtS2[0]);
    dy = (yij - oPtS2[1]);
    dist = dx*dx + dy*dy;
    values[1] = ( std::sqrt(dist) <= rS2) ? 0.1 : 0.0;

    values[2] = 0.0;
    source_->replaceLocalValues(i, values.data());
  }
}

void UnsteadyNonLinAdvDiffReac2dBlockTpetra::setup(){
  createMap();
  initGridAndVel();
  createGraph();
  initFields();
  computeSource();
}

void UnsteadyNonLinAdvDiffReac2dBlockTpetra::assembleMatrix
(const state_type & yState) const
{
   /*
   * Note also that this matrix computes convective and diffusion terms
   * assuming they are on the RHS
   * So: df/dt = -conv + diffusion + else
   */
  auto columnMap = A_->getColMap();
  auto u_ta = u_->getDataNonConst();
  auto v_ta = v_->getDataNonConst();

  int gi{}; int gj{};
  for (LO i=0; i < gptPerProc_; ++i)
  {
    // get the GID of this point
    const auto GID = myGIDs_[i];
    // find global i,j of this point
    globalIDToGiGj(GID, gi, gj);

    auto uij = u_ta[i];
    auto vij = v_ta[i];
    auto A = -2.0 * (eps_ * dxSqInv_ + eps_ * dySqInv_);
    auto B = eps_*dxSqInv_ - uij*dx2Inv_;
    auto F = eps_*dxSqInv_ + uij*dx2Inv_;
    auto D = eps_*dySqInv_ - vij*dy2Inv_;
    auto E = eps_*dySqInv_ + vij*dy2Inv_;

    // Get a view of the current row.
    // Note: I can modify the values, but not the column indices.
    const LO* colInds;
    ST * vals;
    LO numEntries;
    int err = A_->getLocalRowView(i, colInds, vals, numEntries);
    if (err != 0) break;

    for (LO k = 0; k < numEntries; ++k)
    {
      // get the GID of the current column index
      // Note that we need to use the column map to get this
      const auto neighborGID = columnMap->getGlobalElement(colInds[k]);

      // curBlock is a ptr to where current block begins
      ST * const curBlock = &vals[numSpecies_ * numSpecies_ * k];

      if (neighborGID == GID){
	curBlock[0] = A; curBlock[4] = A; curBlock[8] = A;
      }

      // this means the right neighbor
      if (neighborGID == GID+1){
	curBlock[0] = B; curBlock[4] = B; curBlock[8] = B;
      }

      // this means the left neighbor
      if (neighborGID == GID-1){
	curBlock[0] = F; curBlock[4] = F; curBlock[8] = F;
      }

      if (gj>=1 and gj<Ny_-1){
	// this means the neighbor above
	if (neighborGID == GID+Nx_){
	  curBlock[0] = D; curBlock[4] = D; curBlock[8] = D;
	}

	// this means the neighbor above
	if (neighborGID == GID-Nx_){
	  curBlock[0] = E; curBlock[4] = E; curBlock[8] = E;
	}
      }

      if (gj==0 and (neighborGID == GID+Nx_)){
	// check for the neighbor above because of Neumann BC
	// this means the neighbor above
	curBlock[0] = 2.0*eps_*dySqInv_;
	curBlock[4] = 2.0*eps_*dySqInv_;
	curBlock[8] = 2.0*eps_*dySqInv_;
      }

      if (gj==Ny_-1 and (neighborGID == GID-Nx_)) {
	// check for the neighbor below because of Neumann BC
	curBlock[0] = 2.0*eps_*dySqInv_;
	curBlock[4] = 2.0*eps_*dySqInv_;
	curBlock[8] = 2.0*eps_*dySqInv_;
      }
      /* for (LO j = 0; j < numSpecies_; ++j) */
      /* { */
      /* 	const ST newEntry = getCoeff(uij, vij, GID, neighborGID, i); */

      /* 	/\*for (LO i = 0; i < numSpecies_; ++i) {*\/ */
      /* 	LO i = j; */
      /* 	curBlock[i + j * numSpecies_] = newEntry; */
      /* 	/\*}*\/ */
      /* } */
    }
  }//for

  //A_->describe(*out, Teuchos::VERB_EXTREME);
}

void UnsteadyNonLinAdvDiffReac2dBlockTpetra::computeChem( const state_type & C ) const{

  chemReac_->putScalar( static_cast<scalar_type>(0) );
  std::array<ST,3> vals = {};

  // loop over grid points
  for (auto i=0; i<gptPerProc_; i++){
    ST * CC;
    ST * SS;
    // get view of conc entries at this local index
    C.getLocalRowView(i, CC);
    source_->getLocalRowView(i, SS);

    vals[0] = -K_ * CC[0] * CC[1] + SS[0];
    vals[1] = -K_ * CC[0] * CC[1] + SS[1];
    vals[2] = K_*CC[0]*CC[1] - K_*CC[2] + SS[2];

    chemReac_->replaceLocalValues(i, vals.data());
  }
}

}} //namespace rompp::apps
