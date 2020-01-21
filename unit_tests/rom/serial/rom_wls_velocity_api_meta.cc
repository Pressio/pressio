#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
// Currently the WLS System requires calls to ODE package (time discrete residuals)
// and nonlinear solvers (Initial conditions) 
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_WLS"



struct ValidApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;
  using jacobian_type = std::vector<std::vector<scalar_type>>;;
  using dense_matrix_type = std::vector<std::vector<scalar_type>>;

public:
  void velocity(const state_type & y, scalar_type t, velocity_type & f) const
  {};

  velocity_type velocity(const state_type & y, scalar_type t) const{
    velocity_type f;
    return f;
  };

  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const
  {}

  dense_matrix_type applyJacobian(const state_type & y,
				  const dense_matrix_type & B,
				  scalar_type t) const{
    dense_matrix_type A;
    return A;
  }
};

TEST(rom_wls_meta, validVeloAPI){
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( rom::meta::model_meets_velocity_api_for_wls<app_t>::value,"");
}


