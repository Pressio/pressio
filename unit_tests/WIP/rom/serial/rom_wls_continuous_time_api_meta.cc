#include <gtest/gtest.h>
#include "pressio_rom.hpp"

struct ValidApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;
  using jacobian_type = std::vector<std::vector<scalar_type>>;;
  using dense_matrix_type = std::vector<std::vector<scalar_type>>;

public:
  velocity_type createVelocity() const{
    velocity_type f;
    return f;
  };
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    dense_matrix_type A;
    return A;
  }

  void velocity(const state_type & y, scalar_type t, velocity_type & f) const
  {};

  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const
  {}
};

TEST(rom_wls_meta, validVeloAPI){
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( rom::constraints::most_likely_continuous_time_system<app_t>::value,"");
}


