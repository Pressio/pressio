
#include <gtest/gtest.h>
#include "pressio_rom.hpp"
#include "epetra_skeleton.hpp"

TEST(lspg, epetra_types_euler)
{
  using fom_t		= pressio::rom::test::EpetraSkeleton;
  using scalar_t	= typename fom_t::scalar_type;
  using fom_state_t = pressio::containers::Vector<typename fom_t::state_type>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  static_assert(::pressio::rom::concepts::continuous_time_system<fom_t>::value , "");

  using ode_name_t = pressio::ode::implicitmethods::Euler;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_name_t, fom_t, decoder_t, lspg_state_t>::type;
  using lspg_stepper_t = typename lspg_problem::stepper_t;
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
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  static_assert(::pressio::rom::concepts::continuous_time_system<fom_t>::value , "");

  using ode_name_t = pressio::ode::implicitmethods::BDF2;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_name_t, fom_t, decoder_t, lspg_state_t>::type;
  using lspg_stepper_t = typename lspg_problem::stepper_t;
  static_assert(!std::is_void<lspg_stepper_t>::value, "");
}
