
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"
#include "../laplace1dEpetra.hpp"
#include <iostream>
#include <Epetra_MpiComm.h>
#include <Epetra_config.h>
#include <mpi.h>
#include <Epetra_Version.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <stdexcept>
#include <Epetra_CrsMatrix.h>

int main(int argc, char *argv[]){
  using fom_t		= Laplace1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
  using namespace std;
  //---------------------------------------------------------------------------
  // MPI init
  //---------------------------------------------------------------------------
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  //Sets the number of required processors
  assert(Comm.NumProc() == 2);
  //---------------------------------------------------------------------------
  // Initialize Epetra Types
  //---------------------------------------------------------------------------
  shared_ptr<Epetra_Vector> f;
  shared_ptr<Epetra_Vector> u;
  shared_ptr<Epetra_CrsMatrix> A;
  //---------------------------------------------------------------------------
  // Parameters, Setup, and Boundary Conditions
  //---------------------------------------------------------------------------
  vector<double> mu({-1.5, -0.1, 1}); //Parameters: diffusion, advection, expf
  vector<double> domain({0, 2, 0.1}); //1D spatial domain
  vector<double> bc1D({0,0});        //Left and right boundary conditions
  fom_t  appObj(&Comm, mu, domain, bc1D); //Create object
  //Create the linear system
  A = appObj.calculateLinearSystem();  
  f = appObj.calculateForcingTerm();
  u = appObj.createStates();           //states
  //Use AztecOO to solve for the steady states
  appObj.solveForStates(A, u, f);
  //---------------------------------------------------------------------------
  //End and free
  //---------------------------------------------------------------------------
  MPI_Finalize();
  return 0;
}



