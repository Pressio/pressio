#include "CORE_ALL"
#include "APPS_STEADYADVDIFF1D"

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::SteadyAdvDiff1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;

  //---------------------------------------------------------------------------
  // MPI init
  //---------------------------------------------------------------------------
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);
  //---------------------------------------------------------------------------
  // Initialize Epetra Types
  //---------------------------------------------------------------------------
  std::shared_ptr<Epetra_Vector> f;
  std::shared_ptr<Epetra_Vector> u;
  std::shared_ptr<Epetra_CrsMatrix> A;

  //---------------------------------------------------------------------------
  // Parameters, Setup, and Boundary Conditions
  //---------------------------------------------------------------------------
  std::vector<scalar_t> mu({-1, 1, 1}); //Parameters: diffusion, advection, expf
  std::vector<scalar_t> domain({0, 2, 0.05}); //1D spatial domain, xL, xR, dx
  std::vector<scalar_t> bc1D({0,1});        //Left and right boundary conditions
  fom_t  appObj(Comm, mu, domain, bc1D); //Create object

  //---------------------------------------------------------------------------
  //Solve for the states
  //---------------------------------------------------------------------------
  appObj.setup();
  A = appObj.calculateLinearSystem();
  f = appObj.calculateForcingTerm();
  u = appObj.getStates();           //states
  //Use AztecOO to solve for the steady states
  appObj.calculateStates(A, u, f);
  appObj.printStates(u);

  //---------------------------------------------------------------------------
  // Verify Discretization with Manufactured Solution
  //---------------------------------------------------------------------------
  domain[2] = static_cast<scalar_t>(0.0001);
  fom_t  manuFactured(Comm, mu, domain, bc1D); //Create object
  manuFactured.setup();
  A = manuFactured.calculateLinearSystem();
  //Verify that implementation of system is correct
  auto error = manuFactured.verifyImplementation(A);
  assert(error <domain[2]); //Compare to the spatial step size
  std::cout << "L2 Error Manufactured Solution :" << error <<std::endl;

  //---------------------------------------------------------------------------
  //End and free
  //---------------------------------------------------------------------------
  MPI_Finalize();

  return 0;
}
