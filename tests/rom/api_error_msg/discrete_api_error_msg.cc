
#include "pressio_rom.hpp"

struct ValidAppDiscreteTimeApi{
  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using discrete_time_residual_type = state_type;

#if defined FAILMODE5
  using dense_matrix_t = Eigen::MatrixXd;
#else
  using dense_matrix_type = Eigen::MatrixXd;
#endif

public:
    // if FAILMODE1 is defined, exclude this method so to trigger failure
#if defined FAILMODE1
#else
  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
			    pressio::Norm normKind,
			    scalar_type & normR,
  			    Args && ... states) const
  {}
#endif

#if defined FAILMODE2
#else
  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
#if defined FAILMODE5
  				 const dense_matrix_t & B,
  				 dense_matrix_t & A,
#else
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
#endif
  				 Args && ... states) const
  {}
#endif

#if defined FAILMODE3
#else
  discrete_time_residual_type createDiscreteTimeResidual() const
  {
    return discrete_time_residual_type(numDof_);
  }
#endif

#if defined FAILMODE4
#else

#if defined FAILMODE5
  dense_matrix_t createApplyDiscreteTimeJacobianResult(const dense_matrix_t & B) const{
    dense_matrix_t A(numDof_, 3);
    return A;
  }
#else
  dense_matrix_type createApplyDiscreteTimeJacobianResult(const dense_matrix_type & B) const{
    dense_matrix_type A(numDof_, 3);
    return A;
  }
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
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
#if defined FAILMODE5
  using decoder_jac_t	= pressio::containers::MultiVector<typename app_t::dense_matrix_t>;
#else
  using decoder_jac_t	= pressio::containers::MultiVector<typename app_t::dense_matrix_type>;
#endif
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  app_t appobj;
  decoder_jac_t phi(Eigen::MatrixXd(appobj.numDof_,3));
  decoder_t decoderObj(phi);
  typename app_t::state_type yRef;
  lspg_state_t yROM;

  using ode_name_t = pressio::ode::implicitmethods::Arbitrary;
  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_name_t, app_t, lspg_state_t, decoder_t, stepper_order, stepper_n_states>::type;

  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM);
  // std::cout << &lspgProblem << std::endl;

  return 0;
}
