
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "ROM_LSPG_UNSTEADY"

struct ValidApp{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    residual_type & R,
  			    Args & ... states) const
  {
    R.setConstant(1);
    // forward to whatever approriate impl method, e. g.
    // timeDiscreteResidualImpl<step_t>( step, time, f, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  residual_type timeDiscreteResidual(const step_t & step,
  				     const scalar_type & time,
  				     const scalar_type & dt,
  				     Args & ... states) const
  {
    residual_type R(15);
    this->timeDiscreteResidual(step, time, dt, R, std::forward<Args>(states)...);
    return R;
  }

  template <typename step_t, typename ... Args>
  void applyTimeDiscreteJacobian(const step_t & step,
				 const scalar_type & time,
				 const scalar_type & dt,
				 const dense_matrix_type & B,
				 int id,
				 dense_matrix_type & A,
				 Args & ... states) const
  {
    A.setConstant(2);
  }

  template <typename step_t, typename ... Args>
  dense_matrix_type applyTimeDiscreteJacobian(const step_t & step,
  					      const scalar_type & time,
  					      const scalar_type & dt,
  					      const dense_matrix_type & B,
  					      int id,
  					      Args & ... states) const
  {
    dense_matrix_type A(15,3);
    this->applyTimeDiscreteJacobian(step, time, dt, B, id, std::forward<Args>(states)...);
    return A;
  }

};

TEST(rom_lspg, defaultLSPGProblemResidualAPI)
{
  using namespace pressio;
  using app_t    = ValidApp;

  static_assert( rom::meta::model_meets_residual_api_for_unsteady_lspg<app_t>::value,"");
  static_assert( !rom::meta::model_meets_velocity_api_for_unsteady_lspg<app_t>::value,"");

  using scalar_t	= typename app_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<typename app_t::dense_matrix_type>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  app_t appobj;
  decoder_jac_t phi(Eigen::MatrixXd(3,3));
  decoder_t decoderObj(phi);
  typename app_t::state_type yRef;
  lspg_state_t yROM;

  constexpr auto ode_case  = pressio::ode::ImplicitEnum::Arbitrary;

  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

  using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
    pressio::rom::DefaultLSPGUnsteady, ode_case, app_t,
    lspg_state_t, decoder_t, stepper_order, stepper_n_states, scalar_t>;

  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM, 0);
  std::cout << &lspgProblem << std::endl;
}
