
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../../../fom/gold_states_implicit.hpp"

using scalar_t		= double;
using uint_t		= unsigned int;
using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
constexpr auto zero	= ::rompp::utils::constants::zero<scalar_t>();
constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
constexpr auto t0	= zero;

struct FomObserver{
  int snapId_{0};
  eig_dyn_mat A_;
  eig_dyn_vec y0_;

  void resizeRows(int nRows){
    // num of rows = num of dofs
    A_.resize(nRows, 0);
    y0_.resize(nRows,1);
  }

  template <typename T>
  void operator()(size_t step,
		  scalar_t t,
		  const T & y){
    if (step == 0){
      for (auto i=0; i<y0_.rows(); i++)
	y0_[i] = y[i];
    }
    else{
      // snapshots are stored after subtracting init cond
      A_.conservativeResize(Eigen::NoChange,
			    A_.cols()+1);
      for (auto i=0; i<A_.rows(); i++)
	A_(i,snapId_) = y[i]-y0_[i];
      snapId_++;
    }
  }
};


struct FomRunner{
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dEigen;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacobian_t	= typename app_t::jacobian_type;
  using ode_state_t = rompp::containers::Vector<app_state_t>;
  using ode_res_t   = rompp::containers::Vector<app_rhs_t>;
  using ode_jac_t   = rompp::containers::Matrix<app_jacobian_t>;

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

    using stepper_t = rompp::ode::ImplicitStepper<
      ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
    stepper_t stepperObj(y, appobj);

    // define solver
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<
      rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
    rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
    solverO.setTolerance(1e-6);
    solverO.setMaxIterations(100);

    // allocate space in observer to store snapshots
    observer_.resizeRows(totDofs);

    // integrate in time
    rompp::ode::integrateNSteps(stepperObj, y, t0, dt,
				Nsteps, observer_, solverO);

    return y;
  }//run
};


int main(int argc, char *argv[]){
  std::string checkStr {"PASSED"};

  constexpr int Nx = 10, Ny = 10;
  constexpr scalar_t dt = 0.001;
  constexpr auto Nsteps = static_cast<uint_t>(3);

  // run FOM and get snapshots using full mesh case
  // run FOM and collect snapshots
  FomRunner fom(Nx, Ny);
  const auto fomY = fom.run(dt, Nsteps);
  // get snapshots
  const auto & S = fom.getSnapshots();
  //std::cout << S << std::endl;
  // compute SVD to create basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  for (auto i=0; i<U.rows(); ++i){
    if (i % 4 == 0){
      std::cout << std::endl;
      std::cout << i/4 << std::endl;
    }
    std::cout << i << " ";
    for (auto j=0; j<U.cols(); ++j)
      std::cout << std::setprecision(15) << U(i,j) << " ";
    std::cout << std::endl;
  }
  // std::cout << std::setprecision(15) << U << std::endl;
  std::cout << "size(U) " << U.rows() << " " << U.cols() << std::endl;
  std::cout << "Done with snapshots" << std::endl;

  {
    // this is the size of the rom
    constexpr int romSize = Nsteps;

    // the type of the app with sample mesh
    using app_t	 = rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dEigen;

    // create app object
    app_t appobj(Nx, Ny);
    appobj.setup();

    // typedefs used for rom
    using lspg_state_t	= rompp::containers::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::containers::MultiVector<eig_dyn_mat>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // create decoder
    decoder_jac_t phi(U);
    decoder_t decoderObj(phi);

    // my reference state
    const auto yRef = appobj.getInitialState();

    // define ROM state and set to zero
    lspg_state_t yROM(romSize);
    yROM.putScalar(0.0);

    // define LSPG problem
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    // solvers (linear and GN)
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using hessian_t  = rompp::containers::Matrix<eig_dyn_mat>;

    // linear solver
    using solver_tag   = rompp::solvers::linear::iterative::Bicgstab;
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    lin_solver_t linSolverObj;

    // GN solver
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, lin_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-6);
    solver.setMaxIterations(50);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, t0, dt, Nsteps, solver);

    // compute the fom corresponding to our rom final state
    const auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);

    scalar_t rmsErr = zero;
    for (auto i=0; i<yFomFinal.size(); ++i){
      const auto err = yFomFinal[i] - fomY[i];
      rmsErr += err*err;
    }
    auto err = std::sqrt(rmsErr/static_cast<scalar_t>(fomY.size()));
    std::cout << "error = " << err << std::endl;
    if (err > 1e-6)
      checkStr = "FAILED";
  }
  std::cout << checkStr << std::endl;
  return 0;
}
