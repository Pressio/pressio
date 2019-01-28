
#include "CORE_ALL"
#include "ODE_ALL"
#include "SVD_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "QR_BASIC"
#include "../burgers1d_epetra_lspg_euler/lspg_utils.hpp"
#include "../burgers1d_epetra_lspg_euler/burgers1dEpetra.hpp"

const std::vector<double> bdf2Sol =
  {4.9719018997705   ,4.6878137145331   ,3.7881491151657,
   2.3342920452287   , 1.355552338493   ,1.0960796342917,
   1.0595393533521   ,1.0562209827309   ,1.0566764705729,
   1.0579549517173   ,1.0595376412662   ,1.0612369430887,
   1.0635181420314   ,1.0665975221504   ,1.0701653794973,
   1.07356293748   ,1.0764538140891   ,1.0790613676622,
   1.0818312260662   ,1.0850834094093   ,1.0887622944362,
   1.0926539256509   ,1.0965543670099   ,1.1004392803636,
   1.1044023237611   ,1.1085598492751   ,1.1129727267253,
   1.1176261768776   ,1.1224782521064   ,1.1274919395518,
   1.1326674597417   ,1.1380370219587   ,1.1436459848682,
   1.1495149810026   ,1.1556236283773   , 1.161959781159,
   1.1685439733287   ,1.1754010151306   ,1.1825440872759,
   1.1899807566052   ,1.1977183596132   ,1.2057666353442,
   1.2141387640767   ,1.2228508161918   ,1.2319194571235,
   1.2413577532659   ,1.2511762848767   ,1.2613914294993,
   1.27202351768   ,1.2830850813858};


int main(int argc, char *argv[]){

  using fom_t		= Burgers1dEpetra;
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
  assert(Comm.NumProc() == 2);

  //-------------------------------
  // app object
  int numCell = 50; // # fv cells
  std::vector<double> mu({5.0, 0.02, 0.02});
  fom_t appobj(mu, numCell, &Comm);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;

  // store (whichever way you want) the jacobian of the decoder
  constexpr int romSize = 20;
  // store modes computed before from file
  decoder_jac_t phi = rompp::rom::test::getBasis(romSize, numCell, Comm,
  						 appobj.getDataMap());
  const int numBasis = phi.globalNumVectors();
  assert( numBasis == romSize );

  // this is my reference state
  auto & y0n = appobj.getInitialState();
  decoder_t decoderObj(phi);

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  using lspg_problem_types = rompp::rom::DefaultLSPGTypeGenerator<
    fom_t, rompp::ode::ImplicitEnum::BDF2, decoder_t, lspg_state_t>;
  rompp::rom::StepperObjectGenerator<lspg_problem_types> stGen(
      appobj, y0n, decoderObj, yROM, t0);

  using rom_stepper_t = typename lspg_problem_types::rom_stepper_t;

  // GaussNewton solver
  // hessian type: comes up in GN solver, it is (J phi)^T (J phi)
  // Since the rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
  using solver_tag	 = rompp::solvers::linear::LSCG;
  using converged_when_t = rompp::solvers::iterative::default_convergence;
  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
    scalar_t, solver_tag, rompp::solvers::EigenIterative,
    converged_when_t, rom_stepper_t, hessian_t>;
  gnsolver_t solver(stGen.stepperObj_, yROM);
  solver.setTolerance(1e-6);
  solver.setMaxIterations(100);

  // integrate in time
  rompp::ode::integrateNSteps(stGen.stepperObj_, yROM, 0.0, dt, 200, solver);

  // compute the fom corresponding to our rom final state
  using fom_state_w_t = typename lspg_problem_types::fom_state_w_t;
  fom_state_w_t yRf(y0n);
  decoderObj.applyMapping(yROM, yRf);
  yRf += stGen.y0Fom_;
  yRf.data()->Print(std::cout << std::setprecision(14));

  // check against gold solution
  int shift = 0;
  if (rank==1)  shift = 25;
  int myn = yRf.getDataMap().NumMyElements();
  for (auto i=0; i<myn; i++){
    assert(std::abs(yRf[i] - bdf2Sol[i+shift]) < 1e-12 );
   }

  MPI_Finalize();
  return 0;
}
//------------------------------
// end
//------------------------------
