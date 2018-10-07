
#ifndef ODE_REF_APPS_FOR_TESTING_HPP_
#define ODE_REF_APPS_FOR_TESTING_HPP_

namespace rompp{ namespace ode{ namespace testing{ 

struct fakeAppForTraits{
  using scalar_type = double;
  using state_type = std::vector<double>;
  using residual_type = std::vector<double>;

  void residual(const state_type & y,
		residual_type & R,
		scalar_type t) const{
  };
  residual_type residual()const{
    residual_type y;
    return y;
  };
};
//-------------------------------------------
      
struct refAppEigen{

  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;

  void residual(const state_type & y,
		residual_type & R,
		scalar_type t) const{
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++)
      R[i] = y[i];
  };
};

      
}}} // namespace rompp::ode::testing
#endif
