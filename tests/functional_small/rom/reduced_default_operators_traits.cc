
#include <gtest/gtest.h>
#include "pressio/rom_concepts.hpp"

TEST(rom, steady_galerkin)
{
  using namespace pressio;

  using reduced_state_type = Eigen::VectorXd;
  using types = rom::SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using r_t = typename types::reduced_residual_type;
  using j_t = typename types::reduced_jacobian_type;

  static_assert(is_vector_eigen<r_t>::value, "");
  static_assert(is_dense_matrix_eigen<j_t>::value, "");
  static_assert(all_have_traits_and_same_scalar<r_t, j_t>::value, "");
}

TEST(rom, steady_lspg)
{
  using namespace pressio;

  using reduced_state_type = Eigen::VectorXd;
  using types = rom::SteadyLspgDefaultOperatorsTraits<reduced_state_type>;
  using h_t = typename types::hessian_type;
  using g_t = typename types::gradient_type;

  static_assert(is_vector_eigen<g_t>::value, "");
  static_assert(is_dense_matrix_eigen<h_t>::value, "");
  static_assert(all_have_traits_and_same_scalar<h_t, g_t>::value, "");
}
