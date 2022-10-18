
#include <array>

struct MyCustomVector{ 
  std::array<double, 2> d_ = {};
};

struct MyCustomMatrix{ 
  std::array<double, 4> d_ = {};
};

#include "pressio/type_traits.hpp"

namespace pressio{
template<> struct Traits<MyCustomVector>{
  using scalar_type = double;
};

template<> struct Traits<MyCustomMatrix>{
  using scalar_type = double;
};

namespace ops{
MyCustomVector clone(const MyCustomVector & src){ return src; }

MyCustomMatrix clone(const MyCustomMatrix & src){ return src; }

void set_zero(MyCustomVector & o){ o.d_.fill(0); }
void set_zero(MyCustomMatrix & o){ o.d_.fill(0); }

double norm2(const MyCustomVector & v){ 
  return std::sqrt(v.d_[0]*v.d_[0] + v.d_[1]*v.d_[1]); 
}

void update(MyCustomVector & v, double a, const MyCustomVector & v1, double b)
{
  v.d_[0] = v.d_[0]*a + b*v1.d_[0];
  v.d_[1] = v.d_[1]*a + b*v1.d_[1];
}

void scale(MyCustomVector & v, double factor){
  v.d_[0] = v.d_[0]*factor;
  v.d_[1] = v.d_[1]*factor;
}
}}

#include "pressio/solvers_nonlinear.hpp"

struct ValidSystem 
{
  using state_type = MyCustomVector;
  using residual_type = state_type;
  using jacobian_type = MyCustomMatrix;

  state_type createState() const { return state_type{};}
  residual_type createResidual() const { return residual_type{}; }
  jacobian_type createJacobian() const { return jacobian_type{}; }

  void residual(const state_type& x,
                residual_type& res) const
  {
    auto & d = x.d_;
    res.d_[0] =  d[0]*d[0]*d[0] + d[1] - 1.0;
    res.d_[1] = -d[0] + d[1]*d[1]*d[1] + 1.0;   
  }

  void jacobian(const state_type& x, jacobian_type& jac) const 
  {
    auto & J = jac.d_;
    J[0] = 3.0*x.d_[0]*x.d_[0];
    J[1] =  1.0;
    J[2] = -1.0;
    J[3] = 3.0*x.d_[1]*x.d_[1];    
  }
};

struct LinearSolver{
  using matrix_type = MyCustomMatrix;

  void solve(const matrix_type & M, const MyCustomVector & rhs, MyCustomVector & x)
  {
    const auto a = M.d_[0];
    const auto b = M.d_[1];
    const auto c = M.d_[2];
    const auto d = M.d_[3];

    const auto rhs0 = rhs.d_[0];
    const auto rhs1 = rhs.d_[1];

    const auto det = a*d - b*c;
    x.d_[0] = (d*rhs0  - b*rhs1)/det;
    x.d_[1] = (-c*rhs0 + a*rhs1)/det;
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::trace, pressio::log::level::trace});

  using problem_t  = ValidSystem;
  using state_t    = problem_t::state_type;
  problem_t sys;
  state_t y; 
  y.d_[0] = 0.001; 
  y.d_[1] = 0.0001;

  LinearSolver linearSolverObj;
  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(sys, linearSolverObj);
  NonLinSolver.solve(sys, y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y.d_[0] - 1.);
  const auto e2 = std::abs(y.d_[1] - 0.);
  if (e1>1e-8 or e2>1e-8) strOut = "FAILED";

  std::cout <<  strOut << std::endl;
  std::cout << y.d_[0] << " " << y.d_[1] << std::endl;

  pressio::log::finalize();
  return 0;
}
