#include "CONTAINERS_ALL"
#include "APPS_STEADYLINADVDIFF1D"
#include "../fom_gold_states.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){
  int shift = 0;
  if (rank==1) shift = 10;

  for (auto i=0; i<y.MyLength(); i++){
    if (std::abs(y[i] - trueS[i+shift]) > eps) checkStr = "FAILED";
  }
}


int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::SteadyLinAdvDiff1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;

  //--------------------------------------------------
  // MPI init
  //--------------------------------------------------
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);

  //--------------------------------------------------
  // Parameters, Setup, and Boundary Conditions
  //--------------------------------------------------
  //Parameters: diffusion, advection, expf
  std::vector<scalar_t> mu{-1, 1, 1};
  //1D spatial domain, xL, xR
  std::vector<scalar_t> domain{0, 2.0, 0.1};
  //Left and right boundary conditions
  std::vector<scalar_t> bc1D{0, 1.};
  //Create object
  fom_t  appObj(Comm, mu, domain, bc1D);

  //--------------------------------------------------
  //Solve for the states
  //--------------------------------------------------
  appObj.setup();
  appObj.calculateLinearSystem();
  appObj.calculateForcingTerm();
  appObj.solve();
  //  appObj.printState();

  // print LHS matrix
  auto A = appObj.getLHSmatrix();
  A->Print(std::cout << std::setprecision(14));
  std::cout << std::endl;

  // print RHS force vector
  auto f = appObj.getRHSforce();
  f->Print(std::cout << std::setprecision(14));
  std::cout << std::endl;

  auto x = appObj.getGrid();
  auto y = appObj.getState();
  y->Print(std::cout << std::setprecision(14));
  checkSol(rank, *y, pressio::apps::test::steadyAdvDiff1d_N21);

  //--------------------------------------------------
  //End and free
  //--------------------------------------------------
  MPI_Finalize();
  std::cout << checkStr << std::endl;
  return 0;
}
