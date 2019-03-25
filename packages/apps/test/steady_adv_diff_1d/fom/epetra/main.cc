#include "CORE_ALL"
#include "APPS_STEADYADVDIFF1D"

// #include "SOLVERS_NONLINEAR"
// #include "ROM_LSPG"
// #include "QR_BASIC"

/* please make sure you only include things where needed.
 * All the following epetra headers are already included by your class, which
 * is included above via the "APPS_STEADYADVDIFF1D" so it is not needed to
 * include them again. It just creates a burden on the compiler
 */
// #include <Epetra_MpiComm.h>
// #include <Epetra_config.h>
// #include <mpi.h>
// #include <Epetra_Version.h>
// #include <Epetra_Map.h>
// #include <Epetra_Vector.h>
// #include <stdexcept>
// #include <Epetra_CrsMatrix.h>

int main(int argc, char *argv[]){
  using fom_t		= rompp::apps::SteadyAdvDiff1dEpetra;
  using scalar_t	= typename fom_t::scalar_type;
  // using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
  // using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;

  /* be careful with "using namespace".
   * If you use it locally, like here, it is fine,
   * but please NEVER use it in a global scope inside a header file
   * or it will pollute everything. See this for example:
   * https://stackoverflow.com/questions/1452721/why-is-using-namespace-std-considered-bad-practice
 */
  using namespace std;

  /* rather that doing "using namespace", it is preferable to define your own typedefs.
   * Something like:
   * template <typename T> using sh_ptr_t = std::shared_ptr<T>;
   * and
   * using vec_d_t = std::vector<double>;
   *
   * This saves your typing later and clearly tells user what you are doing/using.
   */

  //---------------------------------------------------------------------------
  // MPI init
  //---------------------------------------------------------------------------
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  /* actually, this checks that the # of mpi ranks is = 1.
   * Please change it so that it runs with 2 or 3. */
  //Sets the number of required processors
  assert(Comm.NumProc() == 1);

  //---------------------------------------------------------------------------
  // Initialize Epetra Types
  //---------------------------------------------------------------------------
  shared_ptr<Epetra_Vector> f;
  shared_ptr<Epetra_Vector> u;
  shared_ptr<Epetra_CrsMatrix> Atest;
  shared_ptr<Epetra_CrsMatrix> A;

  //---------------------------------------------------------------------------
  // Parameters, Setup, and Boundary Conditions
  //---------------------------------------------------------------------------
  vector<double> mu({-1, 1, 1}); //Parameters: diffusion, advection, expf
  vector<double> domain({0, 2, 0.05}); //1D spatial domain, xL, xR, dx
  vector<double> bc1D({0,1});        //Left and right boundary conditions
  fom_t  appObj(&Comm, mu, domain, bc1D); //Create object

  //---------------------------------------------------------------------------
  //Solve for the states
  A = appObj.calculateLinearSystem();
  f = appObj.calculateForcingTerm();
  u = appObj.createStates();           //states
  //Use AztecOO to solve for the steady states
  appObj.solveForStates(A, u, f);
  appObj.printStates(u);

  //---------------------------------------------------------------------------
  // Verify Discretization with Manufactured Solution
  //---------------------------------------------------------------------------
  domain[2] = 0.0001;
  fom_t  manuFactured(&Comm, mu, domain, bc1D); //Create object
  Atest = manuFactured.calculateLinearSystem();
  //Verify that implementation of system is correct
  double error = manuFactured.verifyImplementation(Atest);
  std::cout << "L2 Error Manufactured Solution :" << error <<std::endl;

  /* this is still NOT a test. You are just printing to terminal.
   * There is no way to know this test passed because it would require
   * the user to read the output and check visually... which is not feasdible.
   * To test, you should assert some condition for example:
   * assert( error < 1e-12);
   * This will throw an error if that condition is not met.
   */

  //---------------------------------------------------------------------------
  //End and free
  //---------------------------------------------------------------------------
  MPI_Finalize();

  return 0;
}
