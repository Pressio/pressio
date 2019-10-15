
#include <gtest/gtest.h>
#include "ROM_LSPG"
#include "epetra_skeleton.hpp"

TEST(lspg, epetra_types)
{
  using fom_t		= pressio::rom::test::EpetraSkeleton;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  // define LSPG type
  using lspg_problem_types = pressio::rom::DefaultLSPGTypeGenerator<
    fom_t, pressio::ode::ImplicitEnum::Euler, decoder_t, lspg_state_t>;

  using lspg_stepper_t = typename lspg_problem_types::lspg_stepper_t;

  static_assert(!std::is_void<lspg_stepper_t>::value, "");
}
