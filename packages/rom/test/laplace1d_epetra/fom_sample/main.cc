
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"
#include "../laplace1dEpetra.hpp"

int main(int argc, char *argv[]){

  using fom_t		= Laplace1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  // this runs with 3 mpi ranks
  assert(Comm.NumProc() == 3);

  //-------------------------------
  // pseudo steps: 
  // (a) create app object 
  Laplace1dEpetra appObj(&Comm);

  // (b) specify target parameters
  // (c) run problem 
  // (d) this is a test, so we need to check 
  //     that solution is what is expected 

  MPI_Finalize();
  return 0;
}