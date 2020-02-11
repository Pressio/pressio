
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"

namespace{
template <typename result_t>
void readBasis(std::string filename, result_t & phi)
{
  const auto nRows = phi.extent(0);
  const auto nCols = phi.extent(1);

  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, nCols);
  for (int i=0; i<nRows; i++){
    for (int j=0; j<nCols; j++)
      phi(i,j) = A0[i][j];
  }
}
}


template <typename sc_t>
struct myOpsResidualApi
{
  template <typename operand_t>
  static void product(const pressio::apps::arbds::DenseMatrix<sc_t> & A,
  		      const operand_t & b,
		      pressio::apps::arbds::Vector<sc_t> & v)
  {
    // compute: v = A * b
    // b is subscriptable like a regular array, e.g. you can do b[i] or b(i)
    const auto nArows = A.extent(0);
    const auto nAcols = A.extent(1);
    for (auto i=0; i<nArows; ++i)
    {
      v(i) = {};
      for (auto j=0; j<nAcols; ++j)
	v(i) += A(i,j) * b(j);
    }
  }

  static void deep_copy(const pressio::apps::arbds::Vector<sc_t> & from,
			pressio::apps::arbds::Vector<sc_t> & to)
  {
    // here you need do a deep copy from -> to
    to = from;
  }

  static void axpy(sc_t alpha,
		   const pressio::apps::arbds::Vector<sc_t> & x,
		   pressio::apps::arbds::Vector<sc_t> & y)
  {
    // compute y = y + alfa * x
    for (auto i=0; i<y.extent(0); ++i)
      y(i) += alpha * x(i);
  }
};


template <typename dec_jac_t, typename sc_t>
struct myOpsGN
{
  // dot_self overloads needed for the hessian: J^T J
  template <typename result_t>
  static void dot_self(const dec_jac_t & J, result_t & H){
    // J^T J
    for (auto i=0; i<J.extent(1); i++){
      for (auto j=0; j<J.extent(1); j++){
	H(i,j) = pressio::utils::constants::zero<sc_t>();
	for (auto k=0; k<J.extent(0); ++k){
	  H(i,j) += J(k,i) * J(k,j);
	}
      }
    }
  }

  template <typename result_t>
  static result_t dot_self(const dec_jac_t & J){
    result_t result( J.extent(1), J.extent(1) );
    myOpsGN::template dot_self<result_t>(J, result);
    return result;
  }

  // the dot overloads needed for the gradient: J^T R
  template <typename result_t>
  static void dot(const dec_jac_t & J,
		  const pressio::apps::arbds::Vector<sc_t> & r,
		  result_t & g)
  {
    // J^T R
    for (auto i=0; i<J.extent(1); i++){
      g(i) = pressio::utils::constants::zero<sc_t>();
      for (auto j=0; j<r.extent(0); j++){
	g(i) += J(j,i) * r(j);
      }
    }
  }

  static sc_t norm1(const pressio::apps::arbds::Vector<sc_t> & v){
    sc_t result{};
    for (auto i=0; i<v.extent(0); ++i)
      result += std::abs(v(i));
    return result;
  }
  static sc_t norm2(const pressio::apps::arbds::Vector<sc_t> & v){
    sc_t result{};
    for (auto i=0; i<v.extent(0); ++i)
      result += v(i)*v(i);
    return std::sqrt(result);
  }
};


struct EulerLSPGWithResidualApi
{
  using fom_t		= pressio::apps::Burgers1dArbDsResidualApiAdapter;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;

  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  using hessian_t	= pressio::containers::Matrix<eig_dyn_mat>;

  using ops1_t		= myOpsResidualApi<scalar_t>;

  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t,
						      fom_state_t, ops1_t>;

  using opsGN_t = myOpsGN<native_dmat_t, scalar_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithResidualApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    pressio::apps::Burgers1dArbDs appObj(numCell);
    // adapter
    fom_t fomObj(appObj);
    scalar_t dt = 0.01;
    const auto t0 = pressio::utils::constants::zero<scalar_t>();

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes from file
    decoder_jac_t phi(numCell, romSize);
    readBasis("basis.txt", phi);

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    for (auto i=0; i<yRef.extent(0); ++i)
      yRef(i) = pressio::utils::constants::one<scalar_t>();

    // define ROM state
    pressio::containers::ops::resize(yROM_, romSize);
    pressio::containers::ops::fill(yROM_, pressio::utils::constants::zero<scalar_t>());

    // define LSPG type
    using ode_tag = pressio::ode::implicitmethods::Arbitrary;
    using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
    using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

    using lspg_problem	 = pressio::rom::lspg::unsteady::Problem<
      pressio::rom::DefaultLSPGUnsteady, ode_tag, fom_t, lspg_state_t,
      decoder_t, stepper_order, stepper_n_states, scalar_t, ops1_t>;
    using lspg_stepper_t	 = typename lspg_problem::lspg_stepper_t;
    lspg_problem lspgProblem(fomObj, yRef, decoderObj, yROM_, t0);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t, opsGN_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<yFomFinal.extent(0); i++){
      if (std::abs((*yFomFinal.data())(i) - trueY[i]) > 1e-10)
        checkStr = "FAILED";
    }
    std::cout << checkStr <<  std::endl;
  }
};


template <typename nat_dmat_t, typename nat_fom_state_t>
struct myOps2
{
  template <typename operand_t>
  static void product(const nat_dmat_t & A,
  		      const operand_t & vecB,
		      nat_fom_state_t & v)
  {
    // I can do this here bceause I know these are eigen types
    // since this is used below
    v = A * (*vecB.data());
  }
};

struct EulerLSPGWithVelocityApi
{
  using fom_t		= pressio::apps::Burgers1dEigen;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;

  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using hessian_t	= pressio::containers::Matrix<eig_dyn_mat>;

  using ops_t = myOps2<native_dmat_t, native_state_t>;

  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t, ops_t>;

  using ode_tag = pressio::ode::implicitmethods::Euler;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ {};
  
  EulerLSPGWithVelocityApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    Eigen::Vector3d mu(5.0, 0.02, 0.02);
    fom_t appobj( mu, numCell);
    auto t0 = static_cast<scalar_t>(0);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    decoder_jac_t phi =
      pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
    const int numBasis = phi.numVectors();
    if( numBasis != romSize ) throw std::runtime_error("numBasis != romSize");

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    pressio::containers::ops::resize(yROM_, romSize);
    pressio::containers::ops::fill(yROM_, 0.0);

    // define LSPG type
    using lspg_problem	 = pressio::rom::LSPGUnsteadyProblem<
      pressio::rom::DefaultLSPGUnsteady, ode_tag, fom_t, lspg_state_t, decoder_t>;
    using lspg_stepper_t	 = typename lspg_problem::lspg_stepper_t;
    lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM_, t0);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    using gnsolver_t   = pressio::solvers::iterative::GaussNewton<lspg_stepper_t, linear_solver_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM_, linSolverObj);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();

    // this is a reproducing ROM test, so the final reconstructed state
    // has to match the FOM solution obtained with euler, same time-step, for 10 steps
    // const auto trueY = pressio::apps::test::Burg1DtrueImpEulerN20t010;
    const auto trueY = pressio::apps::test::Burgers1dImpGoldStatesBDF1::get(numCell, dt, 0.10);
    for (auto i=0; i<yFomFinal.extent(0); i++){
      if (std::abs(yFomFinal[i] - trueY[i]) > 1e-10)
        checkStr = "FAILED";
    }
    std::cout << checkStr <<  std::endl;
  }
};


int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  EulerLSPGWithVelocityApi LSPGVeloApi;
  const auto veloFomSol = LSPGVeloApi.fomSol_;
  const auto veloRomSol = LSPGVeloApi.yROM_;

  EulerLSPGWithResidualApi LSPGResidApi;
  const auto residFomSol = LSPGResidApi.fomSol_;
  const auto residRomSol = LSPGResidApi.yROM_;

  std::cout << "check that gen coords match" << std::endl;
  // check the reconstructed rom state
  for (auto i=0; i<veloRomSol.extent(0); i++){
    std::cout << std::setprecision(14)
  	      << veloRomSol[i]
  	      << " "
  	      << residRomSol[i]
  	      << std::endl;

    if (std::abs(veloRomSol[i] - residRomSol[i]) > 1e-13)
      checkStr = "FAILED";
  }

  std::cout << "check that fom reconstructed state match" << std::endl;
  // check the reconstructed fom state
  for (auto i=0; i<veloFomSol.size(); i++){
    std::cout << std::setprecision(14)
  	      << veloFomSol[i]
  	      << " "
  	      << residFomSol(i)
  	      << std::endl;

    if (std::abs(veloFomSol[i] - residFomSol(i)) > 1e-13)
      checkStr = "FAILED";
  }

  return 0;
}
