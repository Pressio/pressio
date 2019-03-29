
#include "apps_steady_adv_diff_1d_epetra.hpp"

namespace rompp{ namespace apps{

void SteadyAdvDiff1dEpetra::createMap(){
  //-----------------------------------------------------------------------
  //Create map to be used for the matrices and vectors
  //-----------------------------------------------------------------------
  //Calculate the total number of global nodes
  numGlobalNodes_ = (domain_[1]-domain_[0])/domain_[2] + 1;
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
    if (MyGlobalNodes_[i] ==0 || MyGlobalNodes_[i]== numGlobalNodes_-1)
      NumNz[i] = 1;                             //Boundary conditions
    else
      NumNz[i] = 3;                             //Inner nodes
  //-----------------------------------------------------------------------
  // Create the matrix and vectors to later fill
  //-----------------------------------------------------------------------
  x_ = std::make_shared<nativeVec>(*contigMap_);
  for (int i = 0; i<nodesPerProc_; i++)
    (*x_)[i] = domain_[0]+domain_[2]*MyGlobalNodes_[i];

  //Create forcing term
  f_ = std::make_shared<nativeVec>(*contigMap_);
  //Create A
  A_ = std::make_shared<nativeMatrix>(Copy, *contigMap_, NumNz.get());
  //Create states and leave empty
  u_ = std::make_shared<nativeVec>(*contigMap_);
}


void SteadyAdvDiff1dEpetra::calculateLinearSystem(){
  //-----------------------------------------------------------------------
  //Calculate the components of A
  //-----------------------------------------------------------------------
  //Define individual values
  const int NumEntries = 2;                           //Off diagonal values
  const scalar_type one = static_cast<scalar_type>(1); //Boundary conditions
  const scalar_type two = (-mu_[0]*2.0)/(domain_[2]*domain_[2]); //Main diagonal

  //Define the diagonal values for advection and diffusion
  std::unique_ptr<double[]> Values(new scalar_type[2]);
  std::unique_ptr<int[]> Indices(new int[2]);
  Values[0] = (mu_[0]*1.0)/(domain_[2]*domain_[2])
    -(mu_[1]*1)/(2*domain_[2]);
  Values[1] = (mu_[0]*1.0)/(domain_[2]*domain_[2])
    +(mu_[1]*1)/(2*domain_[2]);

  //-----------------------------------------------------------------------
  //Fill in the values of Matrix A
  //-----------------------------------------------------------------------
  for (int i=0; i<nodesPerProc_; i++){  //Set up the off diagonals
    if (MyGlobalNodes_[i]> 0 && MyGlobalNodes_[i]<numGlobalNodes_-1){
      Indices[0] = MyGlobalNodes_[i] - 1;  //Global Indicies
      Indices[1] = MyGlobalNodes_[i] + 1;  //Global Indicies
      //Fill in the off diagonals
      A_->InsertGlobalValues(MyGlobalNodes_[i], NumEntries, Values.get(),
			     Indices.get());
      //Indicies refer to the column index number here   iagonal 2
      A_->InsertGlobalValues(MyGlobalNodes_[i], 1, &two, MyGlobalNodes_+i);
    }else{
      //Fill in the indicies that refer to the boundary conditions
      A_->InsertGlobalValues(MyGlobalNodes_[i], 1, &one, MyGlobalNodes_+i);
    }
  }
  //Finalize A Matrix Structure
  A_->FillComplete();
}


void SteadyAdvDiff1dEpetra::calculateForcingTerm(){
  /* forcing (3x + x^2) * exp(gamma * x) */
  for (int i = 0; i<nodesPerProc_; i++){
    auto x = (*x_)[i];
    if (MyGlobalNodes_[i]> 0 && MyGlobalNodes_[i]<numGlobalNodes_-1){
      (*f_)[i] = std::exp(mu_[2]*x)*(3.*x + x*x);
    }
    else if (MyGlobalNodes_[i]==0){
      (*f_)[i] = bc1D_[0];
    }
    else if (MyGlobalNodes_[i]==numGlobalNodes_-1){
      (*f_)[i] = bc1D_[1];
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
