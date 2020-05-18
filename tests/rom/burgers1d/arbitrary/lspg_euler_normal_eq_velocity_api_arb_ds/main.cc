
#include "pressio_apps.hpp"
#include "pressio_rom.hpp"
#include "utils_eigen.hpp"

namespace{
template <typename result_t>
void readBasis(std::string filename, result_t & phi)
{
  const auto nRows = phi.extent(0);
  const auto nCols = phi.extent(1);

  std::vector<std::vector<double>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, nCols);
  for (std::size_t i=0; i<nRows; i++){
    for (std::size_t j=0; j<nCols; j++)
      phi(i,j) = A0[i][j];
  }
}
}


template <typename sc_t, typename dec_jac_t>
struct myOps
{
  void time_discrete_residual(sc_t c_R,
			      pressio::apps::arbds::Vector<sc_t> & R,
			      sc_t c_xn,
			      const pressio::apps::arbds::Vector<sc_t> & x_n,
			      sc_t c_xnm1,
			      const pressio::apps::arbds::Vector<sc_t> & x_nm1) const
  {
    // the time discrete residual is: R = c_xn*x_n + c_xnm1*x_nm1 + c_R*f;
    // in input, R contains the velocity so we can do:
    // R = c_xn*x_n + c_xnm1*x_nm1 + c_R*R;
    for (std::size_t i=0; i<R.extent(0); ++i)
      R(i) = c_xn*x_n(i) +c_xnm1*x_nm1(i) + c_R*R(i);
  }

  void time_discrete_jacobian(sc_t c_J,
			      dec_jac_t & J,
			      sc_t c_phi,
			      const dec_jac_t & phi) const
  {
    // this computes the action of the time-disc jac on phi, i.e.: J = tdJac*phi
    // where tdJac = dR/dxn (where R is the time-disc residual above)
    // on input, J contains phi*df/dxn

    for (std::size_t i=0; i<J.extent(0); ++i){
      for (std::size_t j=0; j<J.extent(1); ++j){
    	J(i,j) = c_phi*phi(i,j) + c_J*J(i,j);
      }
    }
  }

  //  y = beta * y + alpha*A*x
  template <typename x_t>
  void product(::pressio::nontranspose mode,
	       const sc_t alpha,
	       const dec_jac_t & A,
	       const x_t & x,
	       const sc_t beta,
	       pressio::apps::arbds::Vector<sc_t> & y) const
  {
    // x is subscriptable like a regular array, e.g. you can do x[i] or x(i)
    const auto nArows = A.extent(0);
    const auto nAcols = A.extent(1);
    for (std::size_t i=0; i<nArows; ++i)
    {
      y(i) = beta*y(i);
      for (std::size_t j=0; j<nAcols; ++j)
	y(i) += alpha * A(i,j) * x(j);
    }
  }

  void deep_copy(pressio::apps::arbds::Vector<sc_t> & to,
		 const pressio::apps::arbds::Vector<sc_t> & from) const
  {
    // here you need do a deep copy from -> to
    to = from;
  }

  void set_zero(pressio::apps::arbds::Vector<sc_t> & vec) const
  {
    for (std::size_t i=0; i<vec.extent(0); ++i)
     vec(i) = static_cast<sc_t>(0);
  }

  void axpy(sc_t alpha,
	    const pressio::apps::arbds::Vector<sc_t> & x,
	    pressio::apps::arbds::Vector<sc_t> & y) const
  {
    // compute y = y + alfa * x
    for (std::size_t i=0; i<y.extent(0); ++i)
      y(i) += alpha * x(i);
  }
};



template <typename dec_jac_t, typename sc_t>
struct myOpsGN
{
  //  C = beta * C + alpha*A^T*B
  template <typename result_t>
  void product(::pressio::transpose modeA,
	       ::pressio::nontranspose modeB,
	       const sc_t alpha,
	       const dec_jac_t & A,
	       const dec_jac_t & B,
	       const sc_t beta,
	       result_t & C) const
  {
    assert( A.extent(0) == B.extent(0) );

    // J^T J
    for (std::size_t i=0; i<A.extent(1); i++){
      for (std::size_t j=0; j<B.extent(1); j++){
    	C(i,j) = beta * C(i,j);
    	for (std::size_t k=0; k<A.extent(0); ++k){
    	  C(i,j) += alpha * A(k,i) * B(k,j);
    	}
      }
    }
  }

  //  C = alpha*A^T*B
  template <typename result_t>
  result_t product(::pressio::transpose modeA,
		   ::pressio::nontranspose modeB,
		   const sc_t alpha,
		   const dec_jac_t & A,
		   const dec_jac_t & B) const
  {
    constexpr auto beta = ::pressio::utils::constants<sc_t>::zero();
    result_t C( A.extent(1), B.extent(1) );
    myOpsGN::product(modeA, modeB, alpha, A, B, beta, C);
    return C;
  }
  // -------------------------------------------------------

  //  y = beta * y + alpha*A^T*x
  template <typename y_t>
  void product(::pressio::transpose mode,
	       const sc_t alpha,
	       const dec_jac_t & A,
	       const pressio::apps::arbds::Vector<sc_t> & x,
	       const sc_t beta,
	       y_t & y) const
  {
    for (std::size_t i=0; i<A.extent(1); i++){
      y(i) = beta * y(i);
      for (std::size_t j=0; j<x.extent(0); j++){
	y(i) += alpha * A(j,i) * x(j);
      }
    }
  }

  // -------------------------------------------------------
  sc_t norm1(const pressio::apps::arbds::Vector<sc_t> & v) const
  {
    sc_t result{};
    for (std::size_t i=0; i<v.extent(0); ++i)
      result += std::abs(v(i));
    return result;
  }

  sc_t norm2(const pressio::apps::arbds::Vector<sc_t> & v) const
  {
    sc_t result{};
    for (std::size_t i=0; i<v.extent(0); ++i)
      result += v(i)*v(i);
    return std::sqrt(result);
  }
};

struct EulerLSPGWithVelocityApi
{
  using fom_t		= pressio::apps::Burgers1dArbDsVelocityApiAdapter;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t   = typename fom_t::dense_matrix_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;

  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  using hessian_t	= pressio::containers::Matrix<eig_dyn_mat>;

  using ops1_t		= myOps<scalar_t, native_dmat_t>;
  using opsGN_t		= myOpsGN<native_dmat_t, scalar_t>;
  static_assert(::pressio::solvers::meta::has_all_needed_methods_for_hessian<
  		opsGN_t, native_dmat_t, hessian_t, scalar_t >::value,
  		"opsGN_t does not have all methods for hessian");

  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t, ops1_t>;

  native_state_t fomSol_ = {};
  lspg_state_t yROM_ = {};

  EulerLSPGWithVelocityApi()
  {
    std::string checkStr {"PASSED"};

    // app object
    constexpr int numCell = 20;
    pressio::apps::Burgers1dArbDs appObj(numCell);
    // adapter
    fom_t fomObj(appObj);
    scalar_t dt = 0.01;
    const auto t0 = pressio::utils::constants<scalar_t>::zero();

    ops1_t myOps1;
    opsGN_t myOps2;

    // read from file the jacobian of the decoder
    constexpr int romSize = 11;
    // store modes from file
    decoder_jac_t phi(numCell, romSize);
    readBasis("basis.txt", phi);

    // create decoder obj
    decoder_t decoderObj(phi, myOps1);

    // for this problem, my reference state = initial state
    native_state_t yRef(numCell);
    for (std::size_t i=0; i<yRef.extent(0); ++i)
      yRef(i) = pressio::utils::constants<scalar_t>::one();

    // define ROM state
    pressio::ops::resize(yROM_, romSize);
    pressio::ops::fill(yROM_, pressio::utils::constants<scalar_t>::zero());

    // define LSPG type
    using ode_tag  = pressio::ode::implicitmethods::Euler;
    using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
      pressio::rom::DefaultLSPGUnsteady, ode_tag, fom_t, lspg_state_t, decoder_t, ops1_t>;
    using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
    lspg_problem lspgProblem(fomObj, yRef, decoderObj, yROM_, t0, myOps1);

    // linear solver
    using solver_tag	 = pressio::solvers::linear::iterative::LSCG;
    using linear_solver_t  = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    // GaussNewton solver
    using gnsolver_t   = pressio::solvers::nonlinear::GaussNewton<lspg_stepper_t, linear_solver_t, opsGN_t>;
    gnsolver_t solver(lspgProblem.getStepperRef(), yROM_, linSolverObj, myOps2);
    solver.setTolerance(1e-13);
    solver.setMaxIterations(4);

    // integrate in time
    pressio::ode::integrateNSteps(lspgProblem.getStepperRef(), yROM_, 0.0, dt, 10, solver);

    // compute the fom corresponding to our rom final state
    auto yFomFinal = lspgProblem.getFomStateReconstructorCRef()(yROM_);
    fomSol_ = *yFomFinal.data();
  }
};

int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  EulerLSPGWithVelocityApi LSPG;
  const auto residFomSol = LSPG.fomSol_;
  const auto residRomSol = LSPG.yROM_;

  std::vector<double> goldFom = {
  				 1.23924053451073,
  				 1.00513782684694,
  				 1.00258750467821,
  				 1.00283530312057,
  				 1.00313333743167,
  				 1.00346287173038,
  				 1.00382706415196,
  				 1.00422955883347,
  				 1.00467438396341,
  				 1.00516599144560,
  				 1.00570930135089,
  				 1.00630975116313,
  				 1.00697335025047,
  				 1.00770673997318,
  				 1.00851726008250,
  				 1.00941302225921,
  				 1.01040299126333,
  				 1.01149707463238,
  				 1.01270622180964,
  				 1.01404253374451
  };

  std::cout << "check that fom reconstructed state match" << std::endl;
  // check the reconstructed fom state
  for (std::size_t i=0; i<residFomSol.extent(0); i++){
    std::cout << std::setprecision(14)
  	      << goldFom[i]
  	      << " "
  	      << residFomSol(i)
  	      << std::endl;

    if (std::abs(goldFom[i] - residFomSol(i)) > 1e-13)
      checkStr = "FAILED";
  }
  std::cout << checkStr <<  std::endl;

  return 0;
}
