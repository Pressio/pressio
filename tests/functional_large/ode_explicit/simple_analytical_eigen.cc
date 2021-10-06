
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include <array>

// constexpr double eps = 1e-12;
// std::string checkStr {"PASSED"};

// template <typename T>
// void checkSol(const T & y, const std::vector<double> & trueS){
//   for (size_t i=0; i< trueS.size(); i++){
//     if (std::abs(y(i) - trueS[i]) > eps) checkStr = "FAILED";
//   }
// }

/*
  let:
    y[0] = a * 0 * t^2
    y[1] = a * 1 * t^2
    y[2] = a * 2 * t^2
    y[3] = a * 3 * t^2
    ...

    dy/dt[0] = 2. * a * 0 * t
    dy/dt[0] = 2. * a * 1 * t
    dy/dt[0] = 2. * a * 2 * t
    dy/dt[0] = 2. * a * 3 * t
    ...

  solve:
     dy/dt = f

  to make this large, we use a state vector of large size
*/

constexpr int num_elements = 50000;

struct MySystem
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using velocity_type = state_type;

  void velocity(const state_type & y,
                scalar_type t,
                velocity_type & rhs) const
  {
    for (int i=0; i<num_elements; ++i)
    {
      const auto my_a = 0.005 * (scalar_type) i;
      rhs(i) = 2. * my_a * t;
    }
  };

  velocity_type createVelocity() const{
    return velocity_type(num_elements);
  };
};

template<class ScalarType>
Eigen::VectorXd gold_solution(ScalarType time)
{
  Eigen::VectorXd g(num_elements);
  for (int i=0; i<num_elements; ++i)
  {
    const auto my_a = 0.005 * (ScalarType) i;
    g(i) = my_a * time*time;
  }
  return g;
}

int main(int argc, char *argv[])
{
#ifndef __INTEL_LLVM_COMPILER

  using app_t		 = MySystem;
  using scalar_t = typename app_t::scalar_type;
  using state_t	 = typename app_t::state_type;

  app_t appObj;
  state_t y(num_elements);

  const scalar_t dt = 0.01;
  auto Nsteps = static_cast<::pressio::ode::step_count_type>(100);

  state_t yssprk3(num_elements);
  {
    yssprk3.setConstant(0.0);
    auto stepperObj = pressio::ode::create_ssprk3_stepper(yssprk3,appObj);
    pressio::ode::advance_n_steps(stepperObj, yssprk3, 0.0, dt, Nsteps);
  }

  state_t yrk4(num_elements);
  {
    yrk4.setConstant(0.0);
    auto stepperObj = pressio::ode::create_rk4_stepper(yrk4,appObj);
    pressio::ode::advance_n_steps(stepperObj, yrk4, 0.0, dt, Nsteps);
  }

  const auto final_time = dt * Nsteps;
  auto gold = gold_solution(final_time);
  scalar_t e1 = {};
  scalar_t e2 = {};
  for (int i=0; i<num_elements; ++i)
  {
    e1 += std::abs( yssprk3(i) - gold(i) );
    e2 += std::abs( yrk4(i) - gold(i) );
  }
  std::cout << std::setprecision(14)
            << std::abs(std::sqrt(e1)) << " "
	    << std::abs(std::sqrt(e2)) << std::endl;

  if ( std::abs(std::sqrt(e1) - 0.0001087430467475) > 1e-12 ){
    std::puts("FAILED");
    return 0;
  }

  if ( std::abs(std::sqrt(e2) - 4.7737329750641e-05) > 1e-12 ){
    std::puts("FAILED");
    return 0;
  }
#endif  

  std::puts("PASSED");
  return 0;
}
