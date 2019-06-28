
#include "ALGEBRA_ALL"
#include "ODE_ALL"
#include "QR_BASIC"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include "../../../fom/gold_states_implicit.hpp"

using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEigen;
using scalar_t		= typename app_t::scalar_type;
using app_state_t	= typename app_t::state_type;
using app_residual_t	= typename app_t::residual_type;
using app_jacobian_t	= typename app_t::jacobian_type;

using ode_state_t = rompp::algebra::Vector<app_state_t>;
using ode_res_t   = rompp::algebra::Vector<app_residual_t>;
using ode_jac_t   = rompp::algebra::Matrix<app_jacobian_t>;

using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using uint_t		= unsigned int;

constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
constexpr auto zero	= ::rompp::algebra::constants::zero<scalar_t>();
constexpr auto t0	= static_cast<scalar_t>(0);

void readMatrixFromFile(std::string filename,
			std::vector<std::vector<double>> & A0,
			int ncols){
  assert( A0.empty() );
  std::ifstream source;
  source.open( filename, std::ios_base::in);
  std::string line, colv;
  std::vector<double> tmpv(ncols);
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (int i=0; i<ncols; i++){
      in >> colv;
      tmpv[i] = atof(colv.c_str());
    }
    A0.emplace_back(tmpv);
  }
  source.close();
}

struct FomObserver{
  int colInd_{0};
  eig_dyn_mat A_;

  void resizeRows(int nRows){
    A_.resize(nRows, 0);
  }

  template <typename T>
  void operator()(size_t step,
		  scalar_t t,
		  const T & y){
    A_.conservativeResize(Eigen::NoChange,
			  A_.cols()+1);
    for (auto i=0; i<A_.rows(); i++)
      A_(i,colInd_) = y[i];
    colInd_++;
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
    solverO.setTolerance(1e-15);
    solverO.setMaxIterations(500);

    // integrate in time
    observer_.resizeRows(totDofs);
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
    using lspg_state_t	= rompp::algebra::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::algebra::MultiVector<eig_dyn_mat>;
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
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_t, ode_case, decoder_t, lspg_state_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0);

    // linear solver
    using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
    using hessian_t  = rompp::algebra::Matrix<eig_dyn_mat>;
    using solver_tag   = rompp::solvers::linear::iterative::Bicgstab;
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    lin_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, lin_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(500);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
    				t0, dt, Nsteps, solver);

    // compute the fom corresponding to our rom final state
    const auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);

    return yFomFinal;
  }//run
};


int main(int argc, char *argv[]){

  /* April 25: this does not work consistenyl
   * for some reason, some runs work, some it does not
   * it might be something going on with eigen? */

  std::string checkStr {"PASSED"};

  // set parameters to use for runs
  constexpr int Nx = 11, Ny = Nx*2-1;
  constexpr scalar_t dt = 0.1;
  constexpr auto Nsteps = static_cast<uint_t>(10);
  constexpr scalar_t fint = Nsteps*dt;

  // // run FOM and collect snapshots
  // FomRunner fom(Nx, Ny);
  // const auto fomY = fom.run(dt, Nsteps);

  // // get snapshots
  // const auto & S = fom.getSnapshots();
  // // do SVD and compute basis
  // Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  // const auto U = svd.matrixU();
  // //std::cout << std::setprecision(15) << U << std::endl;
  // std::cout << S << std::endl;
  // std::cout << "Done with snapshots" << std::endl;

  // std::ofstream file;
  // file.open( "phi.txt" );
  // file << std::fixed
  //      << std::setprecision(15)
  //      << U
  //      << std::endl;
  // file.close();

  constexpr int romSize = Nsteps+1;

  // read gold basis from file
  std::vector<std::vector<double>> U0;
  readMatrixFromFile("basis.txt", U0, romSize);
  Eigen::MatrixXd U(U0.size(), U0[0].size());
  for (size_t i=0; i<U0.size(); i++)
    for (size_t j=0; j<U0[i].size(); j++)
      U(i,j) = U0[i][j];

  // ROM
  LSPGRunner rom(U, Nx, Ny, romSize);
  const auto romY = rom.run(dt, Nsteps);

  {
    // check that solution is right
    using namespace rompp::apps::test;
    const auto trueY
      = NonLinAdvDiffReac2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint);

    scalar_t maxErr={-1.0};
    scalar_t err ={0};
    scalar_t errtmp ={0};
    for (size_t i=0; i<trueY.size(); i++){
      errtmp = (trueY[i] - romY[i]);
      err += errtmp*errtmp;

      if (std::abs(errtmp) > maxErr)
	maxErr = std::abs(errtmp);

      std::cout << std::setprecision(14)
		<< "fom= " << trueY[i]
		<< " romY= " << romY[i]
		<< std::endl;
    }
    const auto rmsErr = std::sqrt(err/(scalar_t) trueY.size());
    std::cout << std::setprecision(14)
	      << " RMSerror= " << rmsErr
	      << " LInferror= " << maxErr << std::endl;
    if ( rmsErr > 1e-6 or maxErr > 1e-6){
      checkStr = "FAILED";
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
