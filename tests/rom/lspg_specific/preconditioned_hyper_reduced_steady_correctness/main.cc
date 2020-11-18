
#include "pressio_rom_lspg.hpp"

struct Preconditioner
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

  void applyPreconditioner(const state_type & yState,
			   residual_type & operand) const
  {
    for (auto i=0; i<operand.size(); ++i) operand(i)+=1.;
  }

  void applyPreconditioner(const state_type & yState,
			   dense_matrix_type & operand) const
  {
    assert(operand.cols()==3);
    for (auto i=0; i<operand.rows(); ++i)
      for (auto j=0; j<operand.cols(); ++j)
	operand(i,j) += static_cast<scalar_type>(1);
  }
};

struct MyFakeApp
{
  int Nst_ = {};
  int Nsm_ = {};
  int Nrom_ = {};
  mutable int countR_ = {};
  mutable int countA_ = {};

public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int Nst, int Nsm, int Nrom)
    : Nst_(Nst), Nsm_(Nsm), Nrom_(Nrom){}

  residual_type createResidual() const{
    state_type v(Nsm_); v.setZero();
    return v;
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    dense_matrix_type A(Nsm_, B.cols()); A.setZero();
    return A;
  }

  void residual(const state_type & state,
		residual_type & f) const
  {
    countR_++;
    if (countR_==1){
      f(0) = 1.1;
      f(1) = 2.2;
      f(2) = 3.3;
      f(3) = 4.4;
      f(4) = 5.5;
    }
    if (countR_==2){
      f(0) = 2.1;
      f(1) = 3.2;
      f(2) = 4.3;
      f(3) = 5.4;
      f(4) = 6.5;
    }
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     dense_matrix_type & A) const
  {
    countA_++;
    if (countA_==1){
      for (auto j=0; j<A.cols(); ++j) A(0,j) = 1.0;
      for (auto j=0; j<A.cols(); ++j) A(1,j) = 2.0;
      for (auto j=0; j<A.cols(); ++j) A(2,j) = 3.0;
      for (auto j=0; j<A.cols(); ++j) A(3,j) = 4.0;
      for (auto j=0; j<A.cols(); ++j) A(4,j) = 5.0;
    }
    if (countA_==2){
      for (auto j=0; j<A.cols(); ++j) A(0,j) = 2.0;
      for (auto j=0; j<A.cols(); ++j) A(1,j) = 3.0;
      for (auto j=0; j<A.cols(); ++j) A(2,j) = 4.0;
      for (auto j=0; j<A.cols(); ++j) A(3,j) = 5.0;
      for (auto j=0; j<A.cols(); ++j) A(4,j) = 6.0;
    }
  }
};

struct MyLinearSolver
{
  using matrix_type = pressio::containers::DenseMatrix<Eigen::MatrixXd>;
  using vec_t       = pressio::containers::Vector<Eigen::VectorXd>;

  std::string & sentinel_;
  int count_ = 0;

  MyLinearSolver(std::string & sentinel) : sentinel_(sentinel){}

  void solve(const matrix_type & A, const vec_t & b, vec_t & x)
  {
    count_++;
    pressio::ops::fill(x, 1.);
    std::cout << "A" << std::endl;
    std::cout << *A.data() << std::endl;
    std::cout << "b" << std::endl;
    std::cout << *b.data() << std::endl;

    if (count_ == 1)
    {
      Eigen::VectorXd bGold(3); bGold.setConstant(-97.);
      if (!bGold.isApprox(*b.data())) sentinel_ = "FAILED";

      Eigen::MatrixXd hGold(3,3); hGold.setConstant(90.);
      if (!hGold.isApprox(*A.data())) sentinel_ = "FAILED";
    }

    if (count_ == 2)
    {
      Eigen::VectorXd bGold(3); bGold.setConstant(-143.5);
      if (!bGold.isApprox(*b.data())) sentinel_ = "FAILED";

      Eigen::MatrixXd hGold(3,3); hGold.setConstant(135.);
      if (!hGold.isApprox(*A.data())) sentinel_ = "FAILED";
    }
  }
};

int main(int argc, char *argv[])
{
  /*
  check that the preconditioned hyp-red lspg problem works as expected.
  We don't solve a real problem but we we use fake data and veryify that
  at every stage involved, from constructor to fom querying to solve,
  the data is supposed to be correct.
  We have the following:
  - stencils mesh: 10 points enumerated as {0,1,...,9}
        o  o  o  o  o  o  o  o  o  o  (points)
        0  1  2  3  4  5  6  7  8  9  (indices of sample mesh points)

  - sample mesh: subset of stencils mesh points {1,4,5,7,8}
        *        *  *     *  *
        0        4  5     7  8  (indices wrt stencil mesh)
        0        1  2     3  4  (enumaration wrt to sample mesh only)

  - romSize = 3
  - start from romState = [0,0,0]
  - use GN with normal equations and do two iterations of the GN solver,
    inner linear solve: simply sets correction to be [1 1 ... 1]
  - linear mapping such that:
              phi = [ 1 1 1;
                      ...  ;
                     10 10 10]

  - the fomObj returns the residual at the sample mesh points:
	f = [1.1 2.2 3.3 4.4 5.5] for first call
	f = [2.1 3.2 4.3 5.4 6.5] for second call

  - the fomObj returns the applyJac which has size (Nsmesh, romSize)
    applyJac = [1 1 1; for first call
                2 2 2;
                3 3 3;
                4 4 4;
                5 5 5]
    applyJac = [2 2 2; for second call
                3 3 3;
                4 4 4;
                5 5 5;
		6 6 6]


  *************************************
  *** first call to linear solver we have ***
  *************************************
  R = applyPrec to f so that we have [2.1 3.2 4.3 5.4 6.5]

  lspgJac after precond =  [2 2 2;
			    3 3 3;
			    4 4 4;
			    5 5 5;
			    6 6 6]

  so that the first call to the linear solver should have:
  b = -lspgJac^T R = [ -97 -97 -97 ]
  neg sign because of the sign convention in pressio

  A = (lspgJac)^T (lspgJac) =
        [90. 90. 90.;
         90. 90. 90.;
         90. 90. 90.]

  *************************************
  *** second call to linear solver we have ***
  *************************************
  R = applyPrec to f so that we have [3.1 4.2 5.3 6.4 7.5]

  lspgJac after precond =  [3 3 3;
			    4 4 4;
			    5 5 5;
			    6 6 6;
			    7 7 7]

  so that the first call to the linear solver should have:
  b = -lspgJac^T R = [ -143.5 -143.5 -143.5 ]
  neg sign because of the sign convention in pressio

  A = (lspgJac)^T (lspgJac) =
        [135. 135. 135.;
         135. 135. 135.;
         135. 135. 135.]
  */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t	= pressio::containers::Vector<native_state_t>;
  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;

  constexpr int Nstencil = 10;
  constexpr int Nsmesh   = 5;
  constexpr int romSize  = 3;

  // app object
  fom_t appObj(Nstencil, Nsmesh, romSize);

  // decoder
  using native_dmatrix_t = typename fom_t::dense_matrix_type;
  using decoder_jac_t	 = pressio::containers::MultiVector<native_dmatrix_t>;
  using decoder_t	 = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
  decoder_jac_t phi(Nstencil, romSize);
  for (auto i=0; i<Nstencil; ++i)
    for (auto j=0; j<romSize; ++j)
      phi(i,j) = (double) (i+1);

  decoder_t  decoderObj(phi);

  // this is my reference state
  native_state_t refState(Nstencil);
  refState.setConstant(1.0);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  // prec obj
  Preconditioner precObj;

  // lspg problem
  auto problem = pressio::rom::lspg::createPreconditionedHyperReducedProblemSteady
    (appObj, decoderObj, romState, refState, precObj);

  // linear solver
  MyLinearSolver linSolverObj(checkStr);

  // GaussNewton solver with normal equations
  auto solver = pressio::rom::lspg::createGaussNewtonSolver(problem, romState, linSolverObj);
  solver.setMaxIterations(2);
  solver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);

  pressio::rom::lspg::solveSteady(problem,romState, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
