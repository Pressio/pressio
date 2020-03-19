
#include "pressio_optimizers.hpp"

struct Rosenbrock1
{
  using scalar_type = double;
  using state_type  = pressio::containers::Vector<Eigen::VectorXd>;

private:
  scalar_type alpha_;

public:
  Rosenbrock1(scalar_type alpha = 100.0) : alpha_(alpha){}

  scalar_type operator()(const state_type & x) const
  {
    const auto n = x.extent(0);
    scalar_type val = {};
    val = alpha_ * std::pow(std::pow( x[0], 2) - x[1], 2) + std::pow( x[0] - 1.0, 2);
    return val;
  }
  void gradient( const state_type & x, state_type &g) const
  {
    g[0] = 2.*alpha_*(std::pow(x[0], 2) - x[1])*2.*x[0] + 2.*(x[0]-1.);
    g[1] = -2.*alpha_*(std::pow(x[0], 2) - x[1]);
  }
};

int main(int argc, char *argv[])
{
  using obj_t   = Rosenbrock1;
  using state_t = typename obj_t::state_type;
  using sc_t    = typename obj_t::scalar_type;

  obj_t rosenbrock;

  // Set Initial Guess
  state_t x(2);
  x[0] = -3.;
  x[1] = -4.;

  using opt_param_t = pressio::optimizers::Parameters<sc_t>;
  opt_param_t MyPars;
  // set parameters

  using opt_prob_t  = pressio::optimizers::Unconstrained<obj_t>;
  opt_prob_t optProblem(MyPars);
  optProblem.solve(rosenbrock, x);

  std::cout << x[0] << " " << x[1] << std::endl;
  return 0;
}
