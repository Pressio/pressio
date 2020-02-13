
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

struct ValidApp{

  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int32_t>;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    residual_type & R,
  			    Args && ... states) const
  {
    R.setConstant(1);
  }

  template <typename step_t, typename ... Args>
  void applyTimeDiscreteJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 int id,
  				 dense_matrix_type & A,
  				 Args && ... states) const
  {
    A.setConstant(2);
  }

  residual_type createTimeDiscreteResidualObject(const state_type & stateIn) const
  {
    residual_type R(numDof_);
    return R;
  }

  dense_matrix_type createApplyTimeDiscreteJacobianObject(const state_type & stateIn,
							  const dense_matrix_type & B) const
  {
    dense_matrix_type A(numDof_, 3);
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
  using native_state_t  = typename app_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t	= pressio::containers::MultiVector<typename app_t::dense_matrix_type>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t>;

  app_t appobj;
  decoder_jac_t phi(Eigen::MatrixXd(appobj.numDof_,3));
  decoder_t decoderObj(phi);
  typename app_t::state_type yRef;
  lspg_state_t yROM;

  using ode_name_t = pressio::ode::implicitmethods::Arbitrary;

  using stepper_order    = ::pressio::ode::types::StepperOrder<1>;
  using stepper_n_states = ::pressio::ode::types::StepperTotalNumberOfStates<2>;

  using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
    pressio::rom::DefaultLSPGUnsteady, ode_name_t, app_t,
    lspg_state_t, decoder_t, stepper_order, stepper_n_states, scalar_t>;

  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM, 0);
  std::cout << &lspgProblem << std::endl;

}
