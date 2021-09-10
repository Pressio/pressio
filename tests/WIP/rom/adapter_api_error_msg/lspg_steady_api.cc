
#include "pressio_rom_lspg.hpp"

struct MyPrec
{
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;
  void applyPreconditioner(const state_type &, residual_type & r) const;
  void applyPreconditioner(const state_type &, dense_matrix_type & jac) const;
};

struct ValidApp{
  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;

public:
#if defined NO_RES
#else
  void residual(const state_type &, residual_type &) const;
#endif

#if defined NO_APP_J
#else
  void applyJacobian(const state_type &,
		     const Eigen::MatrixXd &,
		     Eigen::MatrixXd &) const;
#endif

#if defined NO_RES_C
#else
  residual_type createResidual() const;
#endif

#if defined NO_APP_J_C
#else
  Eigen::MatrixXd createApplyJacobianResult(const Eigen::MatrixXd &) const;
#endif
};

int main(int argc, char *argv[])
{
  using namespace pressio;
  using app_t	      = ValidApp;
  using scalar_t      = typename app_t::scalar_type;
  using native_state_t= typename app_t::state_type;
  using fom_state_t   = pressio::containers::Vector<native_state_t>;

  // these are just randomly set, just for testing
  using lspg_state_t  = pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_jac_t = pressio::containers::MultiVector<Eigen::MatrixXd>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  using lspg_problem =
#if defined PRECOND
    pressio::rom::lspg::impl::composePreconditionedDefaultProblem<
      app_t, decoder_t, lspg_state_t, MyPrec>::type;
#else
    pressio::rom::lspg::impl::composeDefaultProblem<
      app_t, decoder_t, lspg_state_t>::type;
#endif

  // we should never get here because this test fails
  static_assert(std::is_void<lspg_problem>::value, "");

  return 0;
}
