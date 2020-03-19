
#include "pressio_optimizers.hpp"

struct Rosenbrock1
{
  using scalar_type = double;
  using state_type  = pressio::containers::Vector<Eigen::VectorXd>;

private:
  scalar_type alpha_;
  scalar_type const1_;
  scalar_type const2_;

public:
  Rosenbrock1(scalar_type alpha = 100.0)
    : alpha_(alpha), const1_(100.0), const2_(20.0){}

  scalar_type operator()(const state_type & x)
  {
    const auto n = x.extent(0);
    scalar_type val = {};
    for( std::size_t i=0; i<n/2; i++ ) {
      val += alpha_ * std::pow(std::pow( x[2*i], 2) - x[2*i+1], 2);
      val += std::pow( x[2*i] - 1.0, 2);
    }
    return val;
  }
};

int main(int argc, char *argv[])
{
  using obj_t   = Rosenbrock1;
  using state_t = typename obj_t::state_type;
  using sc_t    = typename obj_t::scalar_type;

  obj_t rosenbrock;

  // Set Initial Guess
  state_t x(100);
  for (int i=0; i<50; i++) {
    x[2*i]   = -1.2;
    x[2*i+1] =  1.0;
  }

  using opt_param_t = pressio::optimizers::Parameters<sc_t>;
  opt_param_t MyPars;
  // set parameters

  using opt_prob_t  = pressio::optimizers::Unconstrained<obj_t>;
  opt_prob_t optProblem(MyPars);
  optProblem.solve(rosenbrock, x);

  return 0;
}
