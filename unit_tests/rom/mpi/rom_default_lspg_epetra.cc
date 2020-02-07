
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "ROM_LSPG_UNSTEADY"
#include "epetra_skeleton.hpp"

TEST(lspg, epetra_types_euler)
{
  using fom_t		= pressio::rom::test::EpetraSkeleton;
  using scalar_t	= typename fom_t::scalar_type;
  using fom_state_t = pressio::containers::Vector<typename fom_t::state_type>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t>;

  static_assert(::pressio::rom::meta::model_meets_velocity_api_for_unsteady_lspg<fom_t>::value , "");

  using ode_name_t = pressio::ode::implicitmethods::Euler;
  using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
    pressio::rom::DefaultLSPGUnsteady, ode_name_t, fom_t, lspg_state_t, decoder_t>;
  using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
  static_assert(!std::is_void<lspg_stepper_t>::value, "");
}


TEST(lspg, epetra_types_bdf2)
{
  using fom_t		= pressio::rom::test::EpetraSkeleton;
  using scalar_t	= typename fom_t::scalar_type;
  using fom_state_t = pressio::containers::Vector<typename fom_t::state_type>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, lspg_state_t, fom_state_t>;

  static_assert(::pressio::rom::meta::model_meets_velocity_api_for_unsteady_lspg<fom_t>::value , "");

  using ode_name_t = pressio::ode::implicitmethods::BDF2;
  using lspg_problem = pressio::rom::LSPGUnsteadyProblem<
    pressio::rom::DefaultLSPGUnsteady, ode_name_t, fom_t, lspg_state_t, decoder_t>;
  using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
  static_assert(!std::is_void<lspg_stepper_t>::value, "");
}
