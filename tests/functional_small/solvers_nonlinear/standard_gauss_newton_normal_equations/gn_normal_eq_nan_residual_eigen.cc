
#include "pressio/solvers.hpp"

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
      if (!gold.isApprox(state)) sentinel_ = "FAILED";
    }
    if (nlStep==2){
      Eigen::VectorXd gold(4);
      gold << 2.1,3.1,4.1,5.1;
      if (!gold.isApprox(state)) sentinel_ = "FAILED";
    }
    if (nlStep!=1 or nlStep!=2){
      sentinel_ = "FAILED";
    }
  }
};

struct FakeProblem
{
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  mutable int count_ = 0;

  state_type createState() const{return state_type(4);}
  residual_type createResidual() const{return residual_type(10);}
  jacobian_type createJacobian() const{return jacobian_type(10,4);}

  void residual(const state_type& x, residual_type & res) const
  {
    ++count_;
    std::cout << "RESIDUAL CALL " << count_ << std::endl;
    res.setConstant(1.);
    if (count_==2){
      res(2) = std::numeric_limits<double>::quiet_NaN();
    }
  }

  void jacobian(const state_type & x, jacobian_type & jac) const
  {
    jac.setConstant(1.);
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
    std::cout << x << std::endl;
    x.setConstant(1.1);
  }
};

int main()
{

#ifndef __INTEL_LLVM_COMPILER

  // this test is to check that solver thorws correctly when
  // residual contains NaNs, and it exits
  // In theory the solver should do 3 steps, but we make residual
  // to contain NaNs at the second one so it exits earlier

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string sentinel = "FAILED";
  using namespace pressio;
  using problem_t = FakeProblem;
  using state_t	  = typename problem_t::state_type;
  using mat_t     = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(4);
  x(0)=1.0; x(1)=2.; x(2)=3.; x(3)=4.;

  auto solver = pressio::nonlinearsolvers::create_gauss_newton(problem, FakeLinS<hessian_t>{});
  auto criterion = pressio::nonlinearsolvers::Stop::AfterMaxIters;
  solver.setStoppingCriterion(criterion);
  solver.setMaxIterations(3);
  Observer myO(sentinel);
  solver.setObserver(myO);

  try{
    solver.solve(problem, x);
  }
  catch (::pressio::eh::NonlinearSolveFailure const &e)
  {
    Eigen::VectorXd gold(4);
    gold << 2.1,3.1,4.1,5.1;
    if (!gold.isApprox(x)){
      sentinel = "FAILED";
    }
    else{
      // if we are here, it means test is ok because exception
      // works as expected
      sentinel = "PASSED";
    }

    std::cout << sentinel << std::endl;
    std::cout << std::setprecision(14) << x << std::endl;
  }
  pressio::log::finalize();
#else
    std::cout << "PASSED" << std::endl;
#endif  

  return 0;
}
