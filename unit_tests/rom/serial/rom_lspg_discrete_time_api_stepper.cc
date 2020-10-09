
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

struct ValidApp
{

  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using discrete_time_residual_type = residual_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
  			    Args && ... states) const
  {
    std::cout << "f1\n";
    // R.setConstant(1);
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 dense_matrix_type & A,
  				 Args && ... states) const
  {
    std::cout << "f2\n";
    A.setConstant(2);
    std::cout << "f3\n";
  }

  discrete_time_residual_type createDiscreteTimeResidual() const
  {
    discrete_time_residual_type R(numDof_);
    return R;
  }

  dense_matrix_type createApplyDiscreteTimeJacobianResult(const dense_matrix_type & B) const
  {
    dense_matrix_type A(numDof_, 3);
    return A;
  }
};

TEST(rom_lspg, defaultLSPGProblemResidualAPI)
{
  using namespace pressio;
  using app_t    = ValidApp;
  static_assert( rom::concepts::discrete_time_system<app_t>::value,"");
  static_assert( !rom::concepts::continuous_time_system<app_t>::value,"");

  using scalar_t	= typename app_t::scalar_type;
  using native_state_t  = typename app_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<typename app_t::dense_matrix_type>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  app_t appobj;
  decoder_jac_t phi(Eigen::MatrixXd(appobj.numDof_,3));
  decoder_t decoderObj(phi);
  typename app_t::state_type yRef(appobj.numDof_);
  lspg_state_t yROM(3);

  static_assert(::pressio::rom::concepts::discrete_time_system<app_t>::value, "");

  using ode_name_t = pressio::ode::implicitmethods::Arbitrary;
  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_name_t, app_t, decoder_t, lspg_state_t, stepper_order, stepper_n_states>::type;

  lspg_problem lspgProblem(appobj, decoderObj, yROM, yRef);
  std::cout << &lspgProblem << std::endl;
  // here we just test that the problem is constructed
}
