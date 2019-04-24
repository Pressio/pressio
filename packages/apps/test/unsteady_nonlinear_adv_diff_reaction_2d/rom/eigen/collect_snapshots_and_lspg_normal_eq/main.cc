
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"

using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEigen;
using scalar_t		= typename app_t::scalar_type;
using app_state_t	= typename app_t::state_type;
using app_residual_t	= typename app_t::residual_type;
using app_jacobian_t	= typename app_t::jacobian_type;

using ode_state_t = rompp::core::Vector<app_state_t>;
using ode_res_t   = rompp::core::Vector<app_residual_t>;
using ode_jac_t   = rompp::core::Matrix<app_jacobian_t>;

using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using uint_t		= unsigned int;

constexpr auto zero	= ::rompp::core::constants::zero<scalar_t>();
constexpr auto t0	= static_cast<scalar_t>(0);


struct FomObserver{
  eig_dyn_mat A_;

  void allocate(int nRows, int nCols){
    A_.resize(nRows, nCols);
  }

  template <typename T>
  void operator()(size_t step,
		  double t,
		  const T & y){
    auto j = step;
    for (auto i=0; i<A_.rows(); i++)
      A_(i,j) = y[i];
  }
};


struct FomRunner{
  const int Nx_ = {};
  const int Ny_ = {};
  FomObserver observer_;

  FomRunner(int Nx, int Ny)
    : Nx_{Nx}, Ny_{Ny}, observer_{}{}

  const eig_dyn_mat & getSnapshots() const {
    return observer_.A_;
  }

  ode_state_t run(scalar_t dt, uint_t Nsteps )
  {
    app_t appobj(Nx_, Ny_);
    appobj.setup();
    const auto totDofs = appobj.getUnknownCount();
    const auto y0n = appobj.getInitialState();
    ode_state_t y(y0n);

    constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
    using stepper_t = rompp::ode::ImplicitStepper<
      ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
    stepper_t stepperObj(y, appobj);

    // define solver
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<
      rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
    rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
    solverO.setTolerance(1e-14);
    solverO.setMaxIterations(200);

    // integrate in time
    observer_.allocate(totDofs, Nsteps+1);
    rompp::ode::integrateNSteps(stepperObj, y, t0, dt, Nsteps, observer_, solverO);

    return y;
  }//run
};


struct LSPGRunner{
  const Eigen::MatrixXd & U_;
  const int Nx_ = {};
  const int Ny_ = {};
  const int romSize_ = {};

  LSPGRunner(const Eigen::MatrixXd & U,
	     int Nx, int Ny,
	     int romSize)
    : U_{U}, Nx_{Nx}, Ny_{Ny},
      romSize_{romSize}{}

  ode_state_t run(scalar_t dt, uint_t Nsteps)
  {
    using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::core::MultiVector<eig_dyn_mat>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // app object
    app_t appobj(Nx_, Ny_);
    appobj.setup();

    // create decoder
    decoder_jac_t phi(U_);
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();
    yRef.setConstant(zero);

    // define ROM state and set to zero
    lspg_state_t yROM(romSize_);
    yROM.putScalar(0.0);

    // define LSPG problem
    constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;
    using solver_tag   = rompp::solvers::linear::iterative::LSCG;
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    lin_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<lspg_stepper_t,
								lin_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-14);
    solver.setMaxIterations(200);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
    				t0, dt, Nsteps, solver);

    // compute the fom corresponding to our rom final state
    const auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);

    return yFomFinal;
  }//run
};


constexpr double eps = 1e-10;
std::string checkStr {"PASSED"};

int main(int argc, char *argv[]){

  // set parameters to use for runs
  constexpr int Nx = 11, Ny = Nx*2-1;
  constexpr scalar_t dt = 0.1;
  constexpr scalar_t fint = dt*10;
  constexpr auto Nsteps = static_cast<uint_t>(fint/dt);

  // run FOM and collect snapshots
  FomRunner fom(Nx, Ny);
  const auto fomY = fom.run(dt, Nsteps);

  // get snapshots
  const auto & S = fom.getSnapshots();
  // do SVD and compute basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  std::cout << std::setprecision(7) << S << std::endl;
  std::cout << "Done with snapshots" << std::endl;

  // ROM
  const int romSize = U.cols();
  LSPGRunner rom(U, Nx, Ny, romSize);
  const auto romY = rom.run(dt, Nsteps);

  // this is a reproducing rom, so they should match "exactly"
  for (auto i=0; i<fomY.size(); i++){
    std::cout << "fom= " << fomY[i] << " romY= " << romY[i] << std::endl;
    if (std::abs(fomY[i] - romY[i]) > eps){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
