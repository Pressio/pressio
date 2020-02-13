
#include <gtest/gtest.h>
#include "ROM_BASIC"

// NOTE: here it does not matter to leave all empty since this
// is just for doing type checking
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

TEST(rom_lspg_meta, validVeloAPI){
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( rom::meta::model_meets_velocity_api_for_unsteady_lspg<app_t>::value,"");
  static_assert( rom::meta::is_legitimate_model_for_unsteady_lspg<app_t>::value,"");

  // assert that it does not meet residual api
  static_assert( !rom::meta::model_meets_residual_api_for_unsteady_lspg<app_t>::value,"");
}
