
#include "apps_steady_adv_diff_1d_epetra.hpp"

namespace rompp{ namespace apps{

void SteadyAdvDiff1dEpetra::createMap(){
  //-----------------------------------------------------------------------
  //Create map to be used for the matrices and vectors
  //-----------------------------------------------------------------------
  // total number of dof (we only consider the interior points)
  numGlobalNodes_ = (domain_[1]-domain_[0])/domain_[2] - 1;
  //Make the contiguous Epetra Map with base of 0
  contigMap_ = std::make_shared<Epetra_Map>(numGlobalNodes_, 0, comm_);
  //Double check that the map is actually contiguous
  if (!contigMap_->LinearMap())
    throw std::logic_error("The supposedly contiguous map isn't so.\n");
}

void SteadyAdvDiff1dEpetra::setup(){
  //-----------------------------------------------------------------------
  // Create Map and define indices of nonzero regions for linear system
  //-----------------------------------------------------------------------
  createMap();
  nodesPerProc_= contigMap_->NumMyElements();     //Number of local nodes
  MyGlobalNodes_ = contigMap_->MyGlobalElements();//Retrives the global node
  std::unique_ptr<int[]> NumNz(new int[numGlobalNodes_]);
  for (int i=0; i<nodesPerProc_; i++)
    NumNz[i] = 3;                             //Inner nodes

  //-----------------------------------------------------------------------
  // Create the matrix and vectors to later fill
  //-----------------------------------------------------------------------
  x_ = std::make_shared<nativeVec>(*contigMap_);
  auto dx = domain_[2];
  auto xL = domain_[0];
  for (int i = 0; i<nodesPerProc_; i++)
    (*x_)[i] = xL + dx + dx*MyGlobalNodes_[i];

  //Create forcing term
  f_ = std::make_shared<nativeVec>(*contigMap_);
  //Create A
  A_ = std::make_shared<nativeMatrix>(Copy, *contigMap_, NumNz.get());
  //Create states and leave empty
  u_ = std::make_shared<nativeVec>(*contigMap_);
}


void SteadyAdvDiff1dEpetra::calculateLinearSystem() const{
  //-----------------------------------------------------------------------
  //Calculate the components of A
  //-----------------------------------------------------------------------
  int numEntries{0};
  std::array<scalar_type, 3> Values;
  std::array<int, 3> Indices;
  const scalar_type diagVal = -2.0 * alpha_/(dx_*dx_);
  const scalar_type ip1 = alphaOvDxSq_ + betaOvDx2_;
  const scalar_type im1 = alphaOvDxSq_ - betaOvDx2_;

  for (int i=0; i<nodesPerProc_; i++){
    auto GID = MyGlobalNodes_[i];

    if (GID == 0){
      numEntries = 2;
      Values = {diagVal, ip1, 0.0};
      Indices = {GID, GID+1, 0};
    }
    else if (GID > 0 && GID < numGlobalNodes_-1){
      numEntries = 3;
      Values = {im1, diagVal, ip1};
      Indices = {GID-1, GID, GID+1};
    }
    else /*if (GID == numGlobalNodes_-1)*/{
      numEntries = 2;
      Values = {im1, diagVal, 0.0};
      Indices = {GID-1, GID, 0};
    }

    if (A_->Filled())
      A_->ReplaceGlobalValues(GID, numEntries,
			      Values.data(), Indices.data());
    else
      A_->InsertGlobalValues(GID, numEntries,
			     Values.data(),  Indices.data());
  }

  //Finalize A Matrix Structure
  if(!A_->Filled())
    A_->FillComplete();
}


void SteadyAdvDiff1dEpetra::calculateForcingTerm() const{
  auto uL = bc1D_[0];
  auto uR = bc1D_[1];

  const scalar_type uRf = (alphaOvDxSq_ + betaOvDx2_) * uR;
  const scalar_type uLf = (alphaOvDxSq_ - betaOvDx2_) * uL;

  /* forcing (3x + x^2) * exp(gamma * x) */
  for (int i = 0; i<nodesPerProc_; i++)
  {
    auto GID = MyGlobalNodes_[i];
    auto x = (*x_)[i];
    (*f_)[i] = std::exp(mu_[2]*x)*(3.*x + x*x);

    /* apply correction due to BC */
    if (GID==0){
      (*f_)[i] -= uLf;
    }
    if (GID==numGlobalNodes_-1){
      (*f_)[i] -= uRf;
    }
  }
}

int SteadyAdvDiff1dEpetra::getNumGlobalNodes() const{
  return numGlobalNodes_;
}

std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::getState() const {
  return u_;
}

std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::getGrid() const{
  return x_;
}

void SteadyAdvDiff1dEpetra::solve(){
  //-----------------------------------------------------------------------
  // Set, Create, and Solve Linear System
  //-----------------------------------------------------------------------
  Epetra_LinearProblem Problem(A_.get(), u_.get(), f_.get());
  //Create AztecOO object
  AztecOO Solver(Problem);
  Solver.Iterate(200, 1e-12);
  Solver.NumIters();
  Solver.TrueResidual();
}

void SteadyAdvDiff1dEpetra::printState() const{
  //-----------------------------------------------------------------------
  // Print x and states
  //-----------------------------------------------------------------------
  /* FR: careful with prints, because with MPI things can get messy and
   * prints are not assured to be in order that you think.
   * You can use print methods already implemented in objects like Epetra_Vector
   * where they already made sure prints to work with MPI */

  // for (int i = 0; i<nodesPerProc_; i++){
  //   std::cout<<"x("<<MyGlobalNodes_[i]<<"):\t" << (*x_)[i] <<"\t"
  // 	     <<"u("<<MyGlobalNodes_[i]<<"):\t" << (*u_)[i] <<std::endl;
  // }

  u_->Print( std::cout << std::setprecision(10) );
}

}} //namespace rompp::apps
