

#include "apps_steady_linear_adv_diff_1d_epetra.hpp"

#ifdef HAVE_TRILINOS
namespace pressio{ namespace apps{

void SteadyLinAdvDiff1dEpetra::createMap(){
  //-----------------------------------------------------------------------
  //Create map to be used for the matrices and vectors
  //-----------------------------------------------------------------------
  //Calculate the total number of global nodes
  numGlobalNodes_ = (domain_[1]-domain_[0])/domain_[2] -1 ; //Take out the boundary conditions
  //Make the contiguous Epetra Map with base of 0
  contigMap_ = std::make_shared<Epetra_Map>(numGlobalNodes_, 0, comm_);
  //Double check that the map is actually contiguous
  if (!contigMap_->LinearMap())
    throw std::logic_error("The supposedly contiguous map isn't so.\n");
}

void SteadyLinAdvDiff1dEpetra::setup(){
  //-----------------------------------------------------------------------
  // Create Map and define indices of nonzero regions for linear system
  //-----------------------------------------------------------------------
  createMap();
  nodesPerProc_= contigMap_->NumMyElements();     //Number of local nodes
  MyGlobalNodes_ = contigMap_->MyGlobalElements();//Retrives the global node
  std::unique_ptr<int[]> NumNz(new int[numGlobalNodes_]);
  for (int i=0; i<nodesPerProc_; i++){
    if (MyGlobalNodes_[i] ==0 || MyGlobalNodes_[i]== numGlobalNodes_-1)
      NumNz[i] = 2;                             //Boundary conditions
    else
      NumNz[i] = 3;                             //Inner nodes
  }

  //-----------------------------------------------------------------------
  // Create the matrix and vectors to later fill
  //-----------------------------------------------------------------------
  x_ = std::make_shared<nativeVec>(*contigMap_);
  for (int i = 0; i<nodesPerProc_; i++)
    (*x_)[i] = domain_[2]*(MyGlobalNodes_[i]+1);

  //Create forcing term
  f_ = std::make_shared<nativeVec>(*contigMap_);
  //Create A
  A_ = std::make_shared<nativeMatrix>(Copy, *contigMap_, NumNz.get());
  //Create states and leave empty
  u_ = std::make_shared<nativeVec>(*contigMap_);
}


void SteadyLinAdvDiff1dEpetra::calculateLinearSystem() const{
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
      Values = {{diagVal, ip1, 0.0}};
      Indices = {{GID, GID+1, 0}};
    }
    else if (GID >= 1 && GID <= numGlobalNodes_-2){
      numEntries = 3;
      Values = {{im1, diagVal, ip1}};
      Indices = {{GID-1, GID, GID+1}};
    }
    else /*if (GID == numGlobalNodes_-1)*/{
      numEntries = 2;
      Values = {{im1, diagVal, 0.0}};
      Indices = {{GID-1, GID, 0}};
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


void SteadyLinAdvDiff1dEpetra::calculateForcingTerm() const{
  auto uL = bc1D_[0];
  auto uR = bc1D_[1];
  const scalar_type uRf = (alphaOvDxSq_ + betaOvDx2_) * uR;
  const scalar_type uLf = (alphaOvDxSq_ - betaOvDx2_) * uL;

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

int SteadyLinAdvDiff1dEpetra::getNumGlobalNodes() const{
  return numGlobalNodes_;
}

std::shared_ptr<Epetra_Vector> SteadyLinAdvDiff1dEpetra::getState() const {
  return u_;
}

std::shared_ptr<Epetra_Vector> SteadyLinAdvDiff1dEpetra::getGrid() const{
      return x_;
}

void SteadyLinAdvDiff1dEpetra::solve(){
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

}} //namespace pressio::apps
#endif
