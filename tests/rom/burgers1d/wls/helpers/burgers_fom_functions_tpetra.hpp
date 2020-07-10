
#ifndef PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_TPETRA_HPP_
#define PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_TPETRA_HPP_

#include "utils_tpetra.hpp"

namespace pressio{ namespace testing{ namespace wls{

template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetra & appObj,
		       ::pressio::ode::implicitmethods::Euler odeTag,
		       std::size_t romSize, std::size_t fomSize,
		       rcpcomm_t Comm)
{
  const auto phiNative = pressio::rom::test::tpetra::readBasis("basis_euler.txt", romSize,
							       fomSize, Comm, appObj.getDataMap());
  decoder_d_t decoderObj(phiNative);
  return decoderObj;
}

template <typename decoder_d_t, typename fom_dmat_t, typename rcpcomm_t>
decoder_d_t readBasis( pressio::apps::Burgers1dTpetra & appObj,
		       ::pressio::ode::implicitmethods::BDF2 odeTag,
		       std::size_t romSize, std::size_t fomSize,
		       rcpcomm_t Comm)
{
  const auto phiNative = pressio::rom::test::tpetra::readBasis("basis_bdf2.txt", romSize,
							       fomSize, Comm,
							       appObj.getDataMap());
  decoder_d_t decoderObj(phiNative);
  return decoderObj;
}

template<typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dTpetra & appObj,
		     const y1_t & yFinal, const y2_t & trueY,
		     int rank)
{
  std::string checkStr {"PASSED"};
  auto yFF_v = yFinal.data()->getData();
  int shift = (rank==0) ? 0 : 10;
  const int myn = yFinal.data()->getMap()->getNodeNumElements();
  for (auto i=0; i<myn; i++)
    if (std::abs(yFF_v[i] - trueY[i+shift]) > 1e-10) checkStr = "FAILED";
  return checkStr;
}

}}} //end namespace pressio::testing::wls
#endif
