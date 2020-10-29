
#include "pressio_rom.hpp"

struct MyFakeApp
{
  int Nst_ = {};
  int Nsm_ = {};
  int Nrom_ = {};

public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int Nst, int Nsm, int Nrom)
    : Nst_(Nst), Nsm_(Nsm), Nrom_(Nrom){}

  velocity_type createVelocity() const{
    return state_type(Nsm_);
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    return dense_matrix_type(Nsm_, B.cols());
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    f(0) = 1.1;
    f(1) = 2.2;
    f(2) = 3.3;
    f(3) = 4.4;
    f(4) = 5.5;
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     scalar_type time,
  		     dense_matrix_type & A) const
  {
    for (auto j=0; j<A.cols(); ++j) A(0,j) = 2.0;
    for (auto j=0; j<A.cols(); ++j) A(1,j) = 3.0;
    for (auto j=0; j<A.cols(); ++j) A(2,j) = 4.0;
    for (auto j=0; j<A.cols(); ++j) A(3,j) = 5.0;
    for (auto j=0; j<A.cols(); ++j) A(4,j) = 6.0;
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
      Eigen::VectorXd bGold(3); bGold.setConstant(20.46);
      if (!bGold.isApprox(*b.data())) sentinel_ = "FAILED";

      Eigen::MatrixXd hGold(3,3); hGold.setConstant(158.8);
      if (!hGold.isApprox(*A.data())) sentinel_ = "FAILED";
    }

    if (count_ == 2)
    {
      Eigen::VectorXd bGold(3); bGold.setConstant(-527.34);
      if (!bGold.isApprox(*b.data())) sentinel_ = "FAILED";

      Eigen::MatrixXd hGold(3,3); hGold.setConstant(158.8);
      if (!hGold.isApprox(*A.data())) sentinel_ = "FAILED";
    }
  }
};

int main(int argc, char *argv[])
{
  /*
  this test checks that the hyp-red lspg problem works as expected.
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

  - we do one time steps: t0 -> t1, with dt = 0.2
  - the fomObj returns the velocity at the sample mesh points:
          f = [1.1 2.2 3.3 4.4 5.5]
  - the fomObj returns the applyJac which has size (Nsmesh, romSize)
    applyJac = [2 2 2;
                3 3 3;
                4 4 4;
                5 5 5;
                6 6 6]

  =========================
  =========================
        time step 1
  =========================
  =========================

  *************************************
  *** first call to solver we have ***
  *************************************
  romState     = [0 0 0],
  fomState_n   = [0 ... 0]
  fomState_n-1 = [0 ... 0]

  R = y_n - y_nm-1 - dt*f
    = [-0.22 -0.44 -0.66 -0.88 -1.1]

  lspgJac = I*phi - dt*df/dy*phi
  where I*phi is only taken at the sample mesh points
  and df/dy*phi = [2 2 2;
                   3 3 3;
                   4 4 4;
                   5 5 5;
                   6 6 6]
  we get:
  lspgJac = [2-dt*2 2-dt*2 2-dt*2;    [1.6 1.6 1.6;
             5-dt*3 5-dt*3 5-dt*3;     4.4 4.4 4.4;
             6-dt*4 6-dt*4 6-dt*4; =   5.2 5.2 5.2;
             8-dt*5 8-dt*5 8-dt*5;     7.0 7.0 7.0;
             9-dt*6 9-dt*6 9-dt*6]     7.8 7.8 7.8]

  so that the first call to the linear solver should have:
  b = -lspgJac^T R = [ 20.46 20.46 20.46 ]
  neg sign because of the sign convention in pressio

  A = (lspgJac)^T (lspgJac) =
        [158.8 158.8 158.8;
        158.8 158.8 158.8;
        158.8 158.8 158.8]


  *************************************
  *** second call to solver we have ***
  *************************************
  romState     = [1 1 1],
  fomState_n   = [3 6 9 12 15 18 21 24 27 30]
  fomState_n-1 = [0 ... 0]

  R = y_n - y_nm-1 - dt*f
   = [6-0.22 15-0.44 18-0.66 24-0.88 27-1.1]
   = [5.78 14.56 17.34 23.12 25.9]

  lspgJac = I*phi - dt*df/dy*phi
  where I*phi is only taken at the sample mesh points
  and df/dy*phi = [2 2 2;
                   3 3 3;
                   4 4 4;
                   5 5 5;
                   6 6 6]
  we get:
  lspgJac = [2-dt*2 2-dt*2 2-dt*2;    [1.6 1.6 1.6;
             5-dt*3 5-dt*3 5-dt*3;     4.4 4.4 4.4;
             6-dt*4 6-dt*4 6-dt*4; =   5.2 5.2 5.2;
             8-dt*5 8-dt*5 8-dt*5;     7.0 7.0 7.0;
             9-dt*6 9-dt*6 9-dt*6]     7.8 7.8 7.8]

  so that the first call to the linear solver should have:
  b = -lspgJac^T R = [ -527.34 -527.34 -527.34 ]
  the neg sign because of the sign convention in pressio

  A = (lspgJac)^T (lspgJac) =
  [158.8 158.8 158.8;
   158.8 158.8 158.8;
   158.8 158.8 158.8]
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

  pressio::containers::Vector<Eigen::Matrix<int,-1,1>> hrIndices(Nsmesh);
  hrIndices(0) = 1;
  hrIndices(1) = 4;
  hrIndices(2) = 5;
  hrIndices(3) = 7;
  hrIndices(4) = 8;

  using odetag = pressio::ode::implicitmethods::Euler;
  auto problem = pressio::rom::lspg::createHyperReducedProblemUnsteady<odetag>
    (appObj, decoderObj, romState, refState, hrIndices);

  MyLinearSolver linSolverObj(checkStr);

  // GaussNewton solver with normal equations
  auto solver = pressio::solvers::nonlinear::createGaussNewton
    (problem.stepperRef(), romState, linSolverObj);
  solver.setMaxIterations(2);
  solver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);

  pressio::ode::advanceNSteps(problem.stepperRef(),
			      romState, 0.0, 0.2, 1, solver);

  std::cout << checkStr <<  std::endl;
  return 0;
}
