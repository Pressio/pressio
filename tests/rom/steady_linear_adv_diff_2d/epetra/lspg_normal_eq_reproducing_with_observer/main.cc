
#include "CONTAINERS_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYLINADVDIFF2D"
#include "utils_epetra.hpp"

struct ResidualSampler{
  using vec_t = pressio::containers::Vector<Epetra_Vector>;
  mutable std::vector<double> myR_;

  std::vector<double> getR() const{
    return myR_;
  };

  // void observeResidualEachGNStep(const vec_t & R, int step) const{
  // }

  void observeResidualWhenSolverConverged(const vec_t & R) const{
    auto map = R.data()->Map();
    auto N = map.NumMyElements();
    if (myR_.size() != (size_t) N) myR_.resize(N);
    for (auto i=0; i<N; i++)
      myR_[i] = R[i];
  }
};


int main(int argc, char *argv[]){
  using true_fom_t	= pressio::apps::SteadyLinAdvDiff2dEpetra;
  using fom_adapter_t	= pressio::apps::SteadyLinAdvDiff2dEpetraRomAdapter;
  using scalar_t	= typename fom_adapter_t::scalar_type;
  using native_state	= typename fom_adapter_t::state_type;

  std::string checkStr {"PASSED"};

  // MPI
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 1);

  // we run the FOM and LSPG for same values of parameters
  constexpr scalar_t Pr = 2.48875378592597;
  constexpr scalar_t Re = 41.7029840887032;

  // the discretization to use for solver
  const int Nx = 11, Ny = Nx*2-1;
  // tot num of dof
  const int numDof = (Nx-2)*Ny;

  // -------------------------
  // -------------------------
  // run FOM model first
  // -------------------------
  // -------------------------

  true_fom_t  appObj(Comm, Nx, Ny, Pr, Re);
  appObj.assembleMatrix();
  appObj.fillRhs();
  appObj.solve();
  appObj.printStateToFile("fom.txt");
  pressio::containers::Vector<native_state> yFom(*appObj.getState());

  // -------------------------
  // -------------------------
  // do LSPG ROM
  // -------------------------
  // -------------------------
  using native_state	= typename fom_adapter_t::state_type;
  using fom_state_t = pressio::containers::Vector<native_state>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t>;

  constexpr int romSize = 5;

  // app object for running rom
  fom_adapter_t  appObjROM(Comm, Nx, Ny, Pr, Re);

  // store modes from file
  const decoder_jac_t phi =
    pressio::rom::test::epetra::readBasis("basis.txt", romSize, numDof,
  					 Comm, appObjROM.getDataMap());
  // decoder object
  decoder_t decoderObj(phi);

  // my reference state
  auto yRef = appObjROM.getState();
  yRef->PutScalar( static_cast<scalar_t>(0) );

  // define ROM state and set to zero
  lspg_state_t yROM(romSize);
  yROM.putScalar(0.0);

  // define LSPG type
  using lspg_problem_type = pressio::rom::lspg::steady::DefaultProblemType<
    fom_adapter_t, decoder_t, lspg_state_t>;
  pressio::rom::lspg::steady::ProblemGenerator<lspg_problem_type> lspgProblem(
      appObjROM, *yRef, decoderObj, yROM);
  using rom_system_t = typename lspg_problem_type::lspg_system_t;

  // create sampler for residual (this is passed to solver)
  using observer_t	 = ResidualSampler;
  observer_t myResidSampler;

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::iterative::EigenIterative<
    solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using converged_when_t = pressio::solvers::iterative::default_convergence;
  using gnsolver_t   = pressio::solvers::iterative::GaussNewton<
    rom_system_t, converged_when_t, linear_solver_t, observer_t>;
  gnsolver_t solver(lspgProblem.getSystemRef(), yROM, linSolverObj, myResidSampler);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(200);
  solver.solve(lspgProblem.getSystemRef(), yROM);

  /* the ROM is run for a parameter point that was used to generate
   * the basis, so we should recover the FOM solution exactly */
  auto yFomApprox = lspgProblem.getFomStateReconstructorCRef()(yROM);
  auto errorVec(yFom); 
  pressio::containers::ops::do_update(errorVec, yFom, 1., yFomApprox, -1.);
  const auto norm2err = pressio::containers::ops::norm2(errorVec);
  if( norm2err > 1e-12 ) checkStr = "FAILED";

  // now calculate the residual using the final yROM
  auto fomResidFromFOMState = lspgProblem.getSystemRef().residual(yROM);
  // retrieve the residual that was stored by the sampler
  auto storedR = myResidSampler.getR();
  // the two should match exactly
  for (size_t i=0; i<storedR.size(); i++){
    if( std::abs(storedR[i] - fomResidFromFOMState[i]) > 1e-14)
      checkStr = "FAILED";
  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
