
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

// NOTE: here it does not matter to leave all empty since this
// is just for doing type checking
struct ValidApp
{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using residual_type = state_type;
  using time_discrete_residual_type = residual_type;
  using jacobian_type = std::vector<std::vector<scalar_type>>;
  using dense_matrix_type = std::vector<std::vector<scalar_type>>;

public:

  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
  			    const scalar_type & time,
			    const scalar_type & dt,
  			    time_discrete_residual_type & R,
            pressio::Norm normKind,
            scalar_type & normR,
  			    Args & ... states) const
  {
    // forward to whatever approriate impl method, e. g.
    // timeDiscreteResidualImpl<step_t>( step, time, f, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void applyTimeDiscreteJacobian(const step_t & step,
				 const scalar_type & time,
				 const scalar_type & dt,
				 const dense_matrix_type & B,
				 dense_matrix_type & A,
				 Args & ... states) const
  {}

  time_discrete_residual_type createTimeDiscreteResidual() const
  {
    time_discrete_residual_type R;
    return R;
  }

  dense_matrix_type createApplyTimeDiscreteJacobianResult(const dense_matrix_type & B) const
  {
    dense_matrix_type A;
    return A;
  }

};

TEST(rom_lspg_meta, validResidualAPI){
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( !rom::meta::admissible_system_velocity_api_unsteady_lspg<app_t>::value,"");
  static_assert( rom::meta::admissible_system_time_discrete_residual_api_unsteady_lspg<app_t>::value,"");
}
