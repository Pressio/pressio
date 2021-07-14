
#ifndef PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_KOKKOS_HPP_
#define PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_KOKKOS_HPP_

#include "utils_kokkos.hpp"

namespace pressio { namespace testing { namespace wls {

template <typename decoder_d_t>
decoder_d_t readBasis(pressio::apps::Burgers1dKokkos & appObj,
		      ::pressio::ode::implicitmethods::Euler odeTag,
		      std::size_t romSize, std::size_t fomSize)
{
  using fom_dmat_t = typename decoder_d_t::jacobian_type;
  fom_dmat_t phi("phi", fomSize, 11);
  pressio::rom::test::kokkos::readBasis("basis_euler.txt", romSize, fomSize, *phi.data());
  decoder_d_t decoderObj(std::move(phi));
  return decoderObj;
}

template <typename decoder_d_t>
decoder_d_t readBasis(pressio::apps::Burgers1dKokkos & appObj,
		      ::pressio::ode::implicitmethods::BDF2 odeTag,
		      std::size_t romSize, std::size_t fomSize)
{
  using fom_dmat_t = typename decoder_d_t::jacobian_type;
  fom_dmat_t phi("phi", fomSize, 11);
  pressio::rom::test::kokkos::readBasis("basis_bdf2.txt", romSize, fomSize, *phi.data());
  decoder_d_t decoderObj(std::move(phi));
  return decoderObj;
}

template <typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dKokkos & appObj,
		     const y1_t & yFinal,
		     const y2_t & trueY,
		     std::size_t fomSize)
{
  std::string checkStr{"PASSED"};
  using native_state_t_h = pressio::apps::Burgers1dKokkos::state_type_h;
  native_state_t_h yFinal_h("yFF_h", fomSize);
  Kokkos::deep_copy(yFinal_h, *yFinal.data());

  for(std::size_t i = 0; i < fomSize; i++) {
    if((std::abs(yFinal_h(i) - trueY[i]) > 1e-8) or std::isnan(yFinal_h(i)))
      checkStr = "FAILED";
  }
  return checkStr;
}

}}}//end namespace pressio::testing::wls
#endif
