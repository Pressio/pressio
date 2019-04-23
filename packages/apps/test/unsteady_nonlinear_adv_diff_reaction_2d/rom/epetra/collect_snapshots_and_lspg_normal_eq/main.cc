
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "utils_epetra.hpp"

using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEpetra;
using scalar_t		= typename app_t::scalar_type;
using app_state_t	= typename app_t::state_type;
using app_residual_t	= typename app_t::residual_type;
using ode_state_t	= rompp::core::Vector<app_state_t>;
using ode_res_t		= rompp::core::Vector<app_residual_t>;
using eig_dyn_mat	= Eigen::MatrixXd;
using uint_t		= unsigned int;

constexpr auto zero	= ::rompp::core::constants::zero<scalar_t>();

struct FomObserver{
  int nRows_ = {};
  int nCols_ = {};
  eig_dyn_mat A_;

  FomObserver(int nRows, int nCols)
    : nRows_(nRows),
      nCols_{nCols+1},//because we also have init cond
      A_{nRows_,nCols_}{}

  template <typename T>
  void operator()(size_t step,
  		  double t,
  		  const T & y){
    auto j = step;
    for (auto i=0; i<nRows_; i++) A_(i,j) = y[i];
  }
};


struct FomRunner{
  const Epetra_MpiComm & comm_;
  const int Nx_ = {};
  const int Ny_ = {};

  FomRunner(const Epetra_MpiComm & comm, int Nx, int Ny)
    : comm_{comm}, Nx_{Nx}, Ny_{Ny}{}

  template <typename observer_t>
  ode_state_t run(observer_t & obs, scalar_t dt, uint_t Nsteps ){
    app_t appobj(comm_, Nx_, Ny_);
    appobj.setup();

    const auto y0n = appobj.getInitialState();
    const auto r0n = appobj.residual(y0n, zero);
    ode_state_t y(y0n);
    ode_res_t r(r0n);

    constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
    using stepper_t = rompp::ode::ExplicitStepper<
      ode_case, ode_state_t, app_t, ode_res_t>;
    stepper_t stepperObj(y, appobj, r);

    // integrate in time
    rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs);
    //y.data()->Print(std::cout << std::setprecision(14));
    return y;
  }//run
};


struct LSPGRunner{
  const Eigen::MatrixXd & U_;
  const Epetra_MpiComm & comm_;
  const int Nx_ = {};
  const int Ny_ = {};
  const int romSize_ = {};

  LSPGRunner(const Eigen::MatrixXd & U,
	     const Epetra_MpiComm & comm,
	     int Nx, int Ny,
	     int romSize)
    : U_{U}, comm_{comm},
      Nx_{Nx}, Ny_{Ny},
      romSize_{romSize}{}

  ode_state_t run(scalar_t dt, uint_t Nsteps)
  {
    using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
    using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // app object
    app_t appobj(comm_, Nx_, Ny_);
    appobj.setup();

    // convert basis from eigen to epetra MV
    std::vector<std::vector<scalar_t>> A0;
    std::vector<scalar_t> tmpv(U_.cols());
    for (auto i=0; i<U_.rows(); ++i){
      for (auto j=0; j<U_.cols(); ++j)
	tmpv[j] = U_(i,j);
      A0.emplace_back(tmpv);
    }
    auto phi = ::rompp::apps::test::epetra::convertFromVVecToMultiVec
      (A0, U_.rows(), U_.cols(), comm_, appobj.getDataMap());
    // decoder object
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();
    yRef.PutScalar(zero);

    // define ROM state and set to zero
    lspg_state_t yROM(romSize_);
    yROM.putScalar(0.0);

    auto t0 = static_cast<scalar_t>(0);

    // define LSPG type
    constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;
    using solver_tag   = rompp::solvers::linear::iterative::LSCG;
    using linear_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(200);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
				0.0, dt, Nsteps, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
    //yFomFinal.data()->Print(std::cout << std::setprecision(14));
    return yFomFinal;
  }//run
};



int main(int argc, char *argv[]){

  // fix parameters to use for runs
  constexpr int Nx = 11, Ny = Nx*2-1;
  // dofs = num grid pt * num of species
  constexpr auto numDofs = (Nx-2)*Ny * 3;

  // for time integration
  constexpr scalar_t dt = 0.01;
  constexpr scalar_t fint = dt*10;
  constexpr auto Nsteps = static_cast<uint_t>(fint/dt);

  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 1);

  // run FOM and collect snapshots
  FomObserver obs(numDofs, Nsteps);
  FomRunner fom(Comm, Nx, Ny);
  auto fomY = fom.run(obs, dt, Nsteps);

  // get snapshots
  const auto & S = obs.A_;

  // do SVD and compute basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  std::cout << "Done with snapshots" << std::endl;
  //std::cout << std::setprecision(7) << U << std::endl;

  // for rom
  const int romSize = U.cols();
  LSPGRunner rom(U, Comm, Nx, Ny, romSize);
  auto romY = rom.run(dt, Nsteps);

  for(auto i=0; i<numDofs; ++i)
    std::cout << fomY[i] << " " << romY[i] << std::endl;

  MPI_Finalize();
  return 0;
}











//   // store modes computed before from file
//   const decoder_jac_t phi =
//     rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
//   					 Comm, appObjROM.getDataMap());
//   // decoder object
//   decoder_t decoderObj(phi);

//   // my reference state
//   auto yRef = appObjROM.getState();
//   yRef->PutScalar( static_cast<scalar_t>(0) );

//   // define ROM state and set to zero
//   lspg_state_t yROM(romSize);
//   yROM.putScalar(0.0);

//   // define LSPG type
//   using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
//     fom_adapter_t, decoder_t, lspg_state_t>;
//   rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type> lspgProblem(
//       appObjROM, *yRef, decoderObj, yROM);

//   using rom_system_t = typename lspg_problem_type::lspg_system_t;

//   // linear solver
//   using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
//   using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;
//   using solver_tag   = rompp::solvers::linear::iterative::LSCG;
//   using linear_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
//   linear_solver_t linSolverObj;

//   // GaussNewton solver
//   // hessian comes up in GN solver, it is (J phi)^T (J phi)
//   // rom is solved using eigen, hessian is wrapper of eigen matrix
//   using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
//   using converged_when_t = rompp::solvers::iterative::default_convergence;
//   using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
//     rom_system_t, converged_when_t, linear_solver_t>;
//   gnsolver_t solver(lspgProblem.systemObj_, yROM, linSolverObj);
//   solver.setTolerance(1e-14);
//   solver.setMaxIterations(200);
//   solver.solve(lspgProblem.systemObj_, yROM);

//   /* the ROM is run for a parameter point that was used to generate
//    * the basis, so we should recover the FOM solution exactly */
//   // reconstruct the fom corresponding to our rom final state
//   auto yFomApprox = lspgProblem.yFomReconstructor_(yROM);
//   appObjROM.printStateToFile("rom.txt", *yFomApprox.data());
//   auto errorVec(yFom); errorVec = yFom - yFomApprox;
//   const auto norm2err = rompp::core::ops::norm2(errorVec);
//   if( norm2err > 1e-12 ) checkStr = "FAILED";

//   std::cout << std::setprecision(15) << norm2err << std::endl;

//   MPI_Finalize();
//   std::cout << checkStr <<  std::endl;
//   return 0;
// }
