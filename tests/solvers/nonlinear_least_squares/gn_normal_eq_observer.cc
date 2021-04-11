
#include "pressio_solvers.hpp"

struct Observer
{
  std::string & sentinel_;

  Observer(std::string & sentinel) : sentinel_(sentinel){}

  template<typename state_t>
  void operator()(int nlStep, const state_t & state) const
  {
    PRESSIOLOG_DEBUG("NONLIN observer");
    std::cout << nlStep << " " << std::endl;
    std::cout << *state.data() << std::endl;

    if (nlStep==1){
      Eigen::VectorXd gold(4);
      gold << 1.,2.,3.,4.;
      if (!gold.isApprox(*state.data())) sentinel_ = "FAILED";
    }
    if (nlStep==2){
      Eigen::VectorXd gold(4);
      gold << 2.1,3.1,4.1,5.1;
      if (!gold.isApprox(*state.data())) sentinel_ = "FAILED";
    }
    if (nlStep==3){
      Eigen::VectorXd gold(4);
      gold << 3.2,4.2,5.2,6.2;
      if (!gold.isApprox(*state.data())) sentinel_ = "FAILED";
    }
  }
};

struct FakeProblem
{
  using scalar_type	= double;
  using state_type	= pressio::containers::Vector<Eigen::VectorXd>;
  using residual_type	= state_type;
  using jacobian_type	= pressio::containers::DenseMatrix<Eigen::MatrixXd>;

  residual_type createResidual() const{return residual_type(10);}
  jacobian_type createJacobian() const{return jacobian_type(10,4);}

  void residual(const state_type& x, residual_type & res) const
  {
    res.data()->setConstant(1.);
  }

  void jacobian(const state_type & x, jacobian_type & jac) const
  {
    jac.data()->setConstant(1.);
  }
};

template<class T>
struct FakeLinS
{
  int count_=0;
  using matrix_type = T;

  template<typename A_t, typename b_t, typename x_t>
  void solve(const A_t & A, const b_t & b, x_t & x)
  {
    ++count_;
    std::cout << *x.data() << std::endl;
    x.data()->setConstant(1.1);
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string sentinel = "PASSED";
  using namespace pressio;
  using problem_t = FakeProblem;
  using state_t	  = typename problem_t::state_type;
  using mat_t     = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(4);
  x(0)=1.0; x(1)=2.; x(2)=3.; x(3)=4.;

  auto solver = pressio::solvers::nonlinear::createGaussNewton(problem, x, FakeLinS<hessian_t>{});
  auto criterion = pressio::solvers::nonlinear::stop::afterMaxIters;
  solver.setStoppingCriterion(criterion);
  solver.setMaxIterations(3);
  Observer myO(sentinel);
  solver.setObserver(myO);
  solver.solve(problem, x);

  Eigen::VectorXd gold(4);
  gold << 3.2,4.2,5.2,6.2;
  if (!gold.isApprox(*x.data())) sentinel = "FAILED";

  std::cout << sentinel << std::endl;
  std::cout << std::setprecision(14) << *x.data() << std::endl;
}
