
#include "pressio_rom_galerkin.hpp"

struct ValidApp{
  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using velocity_type	= state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
#if defined NO_VEL
#else
  void velocity(const state_type &,
		const scalar_type &,
		velocity_type &) const;
#endif

#if defined NO_VEL_C
#else
  velocity_type createVelocity() const;
#endif
};

int main(int argc, char *argv[])
{
  using namespace pressio;
  using app_t	      = ValidApp;
  using scalar_t      = typename app_t::scalar_type;
  using native_state_t= typename app_t::state_type;
  using native_dense_mat = typename app_t::dense_matrix_type;
  using fom_state_t   = pressio::containers::Vector<native_state_t>;

  // these are just randomly set, just for testing
  using lspg_state_t  = pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_jac_t = pressio::containers::MultiVector<native_dense_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  using ode_name_t = pressio::ode::explicitmethods::Euler;
  using problem_t = pressio::rom::galerkin::composeDefaultProblem_t<
    ode_name_t, app_t, decoder_t, lspg_state_t>;

  // we should never get here because this test fails
  static_assert(std::is_void<problem_t>::value, "");

  return 0;
}
