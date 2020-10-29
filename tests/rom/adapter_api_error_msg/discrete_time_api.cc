
#include "pressio_rom_galerkin.hpp"
#include "pressio_rom_lspg.hpp"

struct ValidAppDiscreteTimeApi{

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using discrete_time_residual_type = state_type;

#if defined NO_MAT_ALIAS
  using dense_matrix_t = Eigen::MatrixXd;
#else
  using dense_matrix_type = Eigen::MatrixXd;
#endif

public:
#if defined NO_RES
#else
  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
  			    Args && ... states) const;
#endif

#if defined NO_APP_J
#else
  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
#if defined NO_MAT_ALIAS
  				 const dense_matrix_t & B,
  				 dense_matrix_t & A,
#else
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
#endif
  				 Args && ... states) const;
#endif

#if defined NO_RES_C
#else
  discrete_time_residual_type createDiscreteTimeResidual() const;
#endif

#if defined NO_APP_J_C
#else

#if defined NO_MAT_ALIAS
  dense_matrix_t createApplyDiscreteTimeJacobianResult(const dense_matrix_t & B) const;
#else
  dense_matrix_type createApplyDiscreteTimeJacobianResult(const dense_matrix_type & B) const;
#endif

#endif
};

int main(int argc, char *argv[])
{
  using namespace pressio;
  using app_t    = ValidAppDiscreteTimeApi;
  using scalar_t	= typename app_t::scalar_type;
  using native_state_t  = typename app_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using native_mat_t =
#if defined NO_MAT_ALIAS
    typename app_t::dense_matrix_t;
#else
    typename app_t::dense_matrix_type;
#endif

  using decoder_jac_t	= pressio::containers::MultiVector<native_mat_t>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  using ode_name_t = pressio::ode::implicitmethods::Arbitrary;
  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;


  using problem =
#if defined DO_LSPG
    pressio::rom::lspg::composeDefaultProblem<
#else
    pressio::rom::galerkin::composeDefaultProblem<
#endif
    ode_name_t, app_t, decoder_t, rom_state_t, stepper_order, stepper_n_states>::type;

  // we should never get here because this test fails
  static_assert(std::is_void<problem>::value, "");

  return 0;
}
