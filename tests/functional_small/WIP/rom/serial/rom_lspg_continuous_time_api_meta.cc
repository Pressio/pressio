
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

struct ValidApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;
  using jacobian_type = std::vector<std::vector<scalar_type>>;
  using dmat_type = std::vector<std::vector<scalar_type>>;

public:
  void velocity(const state_type & y, scalar_type t, velocity_type & f) const
  {};

  void applyJacobian(const state_type & y,
         const dmat_type & B,
         scalar_type t,
         dmat_type & A) const
  {}

  velocity_type createVelocity() const{
    velocity_type f;
    return f;
  };

  dmat_type createApplyJacobianResult(const dmat_type & B) const{
    dmat_type A;
    return A;
  }
};

TEST(rom_lspg_meta, validVeloAPI)
{
  using app_t    = ValidApp;
  using namespace pressio;
  using dmat_type = std::vector<std::vector<double>>;

  static_assert( rom::constraints::continuous_time_system_with_user_provided_apply_jacobian
    <app_t, dmat_type>::value,"");
  static_assert( !rom::constraints::discrete_time_system_with_user_provided_apply_jacobian
    <app_t, dmat_type>::value,"");
}
