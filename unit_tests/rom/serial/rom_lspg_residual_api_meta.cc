
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ROM_LSPG"

// NOTE: here it does not matter to leave all empty since this
// is just for doing type checking
struct ValidApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using residual_type = state_type;
  using dense_matrix_type = std::vector<std::vector<scalar_type>>;

private:
  template <typename step_t>
  void timeDiscreteResidualImpl(const step_t & step,
				const scalar_type & time,
				residual_type & f,
				const state_type & y1,
				const state_type & y2) const
  {}

public:

  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
  			    const scalar_type & time,
  			    residual_type & f,
  			    Args & ... states) const
  {
    // forward to whatever approriate impl method, e. g.
    // timeDiscreteResidualImpl<step_t>( step, time, f, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  residual_type timeDiscreteResidual(const step_t & step,
  				     const scalar_type & time,
  				     Args & ... states) const
  {
    residual_type g;
    this->timeDiscreteResidual(step, time, g, std::forward<Args>(states)...);
    return g;
  }

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

TEST(rom_lspg_meta, validResidualAPI){
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( rom::meta::model_meets_residual_api_for_lspg<app_t>::value,"");
  static_assert( rom::meta::is_legitimate_model_for_lspg<app_t>::value,"");

  // assert that it does not meet velocity api
  static_assert( !rom::meta::model_meets_velocity_api_for_lspg<app_t>::value,"");
}
