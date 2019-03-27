#include "CORE_ALL"
#include "APPS_STEADYADVDIFF1D"
#include "../fom_gold_states.hpp"

constexpr double eps = 1e-12;

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS){
  if (rank==0){
    assert(std::abs(y[0] - trueS[0]) < eps);
    assert(std::abs(y[1] - trueS[1]) < eps);
    assert(std::abs(y[2] - trueS[2]) < eps);
    assert(std::abs(y[3] - trueS[3]) < eps);
    assert(std::abs(y[4] - trueS[4]) < eps);
    assert(std::abs(y[5] - trueS[5]) < eps);
    assert(std::abs(y[6] - trueS[6]) < eps);
    assert(std::abs(y[7] - trueS[7]) < eps);
    assert(std::abs(y[8] - trueS[8]) < eps);
    assert(std::abs(y[9] - trueS[9]) < eps);
    assert(std::abs(y[10] - trueS[10]) < eps);
  }
  if (rank==1){
    assert(std::abs(y[0] - trueS[11]) < eps);
    assert(std::abs(y[1] - trueS[12]) < eps);
    assert(std::abs(y[2] - trueS[13]) < eps);
    assert(std::abs(y[3] - trueS[14]) < eps);
    assert(std::abs(y[4] - trueS[15]) < eps);
    assert(std::abs(y[5] - trueS[16]) < eps);
    assert(std::abs(y[6] - trueS[17]) < eps);
    assert(std::abs(y[7] - trueS[18]) < eps);
    assert(std::abs(y[8] - trueS[19]) < eps);
    assert(std::abs(y[9] - trueS[20]) < eps);
  }
}


int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::SteadyAdvDiff1dEpetra;
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
  std::vector<scalar_t> bc1D{0,1};
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

  auto y = appObj.getState();
  auto x = appObj.getGrid();
  y->Print(std::cout << std::setprecision(14));
  checkSol(rank, *y, rompp::apps::test::steadyAdvDiff_N20);

  //--------------------------------------------------
  //End and free
  //--------------------------------------------------
  MPI_Finalize();

  return 0;
}
