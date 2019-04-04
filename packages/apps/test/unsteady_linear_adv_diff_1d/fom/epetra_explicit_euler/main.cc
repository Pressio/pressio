
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYLINADVDIFF1D"

int main(int argc, char *argv[]){
  
  using app_t         =rompp::apps::UnsteadyLinAdvDiff1dEpetra;
  using scalar_t      =typename app_t::scalar_type;
  using app_state_t   =typename app_t::state_type;
  using app_residual_t = typename app_t::residual_type;  //Need residual for unsteady case?

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 2);
  
  std::vector<scalar_t> mu{-1, 1, 1};  
  std::vector<scalar_t> domain{0, 2.0, 0.1};
  std::vector<scalar_t> bc1D{0, 1.0};
 
  app_t appObj(Comm, mu, domain, bc1D);
  appObj.setup();
  auto y0n = appObj.getState();  //Get the initial states
  
  // auto r0n = appObj.residual(y0n, r0n);
  // std::cout << "HELLO WORLD" << std::endl;

  MPI_Finalize();
  return 0;
}
