
#include "pressio_rom_galerkin.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"


struct GalerkinBDF1WithResidualApi
{
  using fom_t		= pressio::apps::Burgers1dEigenDiscreteTimeApi;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat   = Eigen::Matrix<scalar_t, -1, -1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using rom_jacobian_t = pressio::containers::DenseMatrix<eig_dyn_mat>;

  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  native_state_t fomSol_ = {};
  rom_state_t yROM_ = {};

  GalerkinBDF1WithResidualApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    Eigen::Vector3d mu(5.0, 0.02, 0.02);
    fom_t appobj( mu, numCell);
    scalar_t dt = 0.01;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes computed before from file
    decoder_jac_t phi =
      pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell);
    const int numBasis = phi.numVectors();
    if( numBasis != romSize )
      throw std::runtime_error("numBasis != romSize");

    // create decoder obj
    decoder_t decoderObj(phi);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    yRef.setConstant(1);

    // define ROM state
    ::pressio::ops::resize(yROM_, romSize);
    ::pressio::ops::fill(yROM_, 0.0);

    //using ode_tag = pressio::ode::implicitmethods::Arbitrary;
    // using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
    // using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
    // using problem_t = pressio::rom::galerkin::composeDefaultProblem<
    //   ode_tag, fom_t, decoder_t, rom_state_t, rom_jacobian_t,
    // stepper_order, stepper_n_states>::type;
    // problem_t Problem(appobj, decoderObj, yROM_, yRef);
    auto Problem =
      pressio::rom::galerkin::createDefaultProblem<rom_jacobian_t, 1, 2>
      (appobj, decoderObj, yROM_, yRef);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, rom_jacobian_t>;
    linear_solver_t linSolverObj;

    // nonlinear system
    auto solver = pressio::rom::galerkin::createNewtonRaphsonSolver(Problem, yROM_, linSolverObj);
    solver.setTolerance(1e-12);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::rom::galerkin::solveNSteps(Problem, yROM_, 0.0, dt, 15, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = Problem.fomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};

int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  GalerkinBDF1WithResidualApi GalerkinResidApi;
  const auto residFomSol = GalerkinResidApi.fomSol_;
  const auto residRomSol = GalerkinResidApi.yROM_;

  std::cout << "check that fom reconstructed state match" << std::endl;

  const std::vector<double> gold = {
1.35670784589790e+00,
1.00980511428193e+00,
1.00391615986692e+00,
1.00425132885340e+00,
1.00469776316420e+00,
1.00519182353963e+00,
1.00573784922521e+00,
1.00634130030401e+00,
1.00700821600573e+00,
1.00774527082009e+00,
1.00855984111213e+00,
1.00946007894552e+00,
1.01045499373951e+00,
1.01155454237401e+00,
1.01276972872420e+00,
1.01411271389327e+00,
1.01559693785592e+00,
1.01723725390909e+00,
1.01905007725692e+00,
1.02105354928633e+00};

  if ((std::size_t)residFomSol.size() != gold.size())
    checkStr = "FAILED";

  // check the reconstructed fom state
  for (auto i=0; i<residFomSol.size(); i++){
    std::cout << std::setprecision(14)
  	      << gold[i]
  	      << " "
  	      << residFomSol[i]
  	      << std::endl;

    if (std::abs(gold[i] - residFomSol[i]) > 1e-11){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr <<  std::endl;
  return 0;
}
