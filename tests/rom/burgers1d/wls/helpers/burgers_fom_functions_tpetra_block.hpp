
#ifndef PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_TPETRA_BLOCK_HPP_
#define PRESSIO_TESTS_WLS_BURGERS1D_FOM_FUNCTIONS_TPETRA_BLOCK_HPP_

#include "utils_tpetra.hpp"

namespace pressio { namespace testing { namespace wls {

template <typename decoder_d_t, typename rcpcomm_t>
decoder_d_t readBasis(pressio::apps::Burgers1dTpetraBlock & appObj,
		      ::pressio::ode::implicitmethods::Euler odeTag,
		      std::size_t romSize, std::size_t fomSize,
		      rcpcomm_t Comm)
{
  // read basis in a tpetra mv
  auto tp_phi = pressio::rom::test::tpetra::readBasis("basis_euler.txt", romSize, fomSize,
						      Comm, appObj.getDataMap());

  // construct a tpetra_block mv from a tpetra mv
  using tpb_mv = Tpetra::BlockMultiVector<>;
  tpb_mv tpb_phi(*tp_phi.data(), *appObj.getDataMap(), 1 /*1 is block size*/);
  typename decoder_d_t::jacobian_type phi(tpb_phi);

  decoder_d_t decoderObj(phi);
  return decoderObj;
}

template <typename decoder_d_t, typename rcpcomm_t>
decoder_d_t readBasis(pressio::apps::Burgers1dTpetraBlock & appObj,
		      ::pressio::ode::implicitmethods::BDF2 odeTag,
		      std::size_t romSize, std::size_t fomSize,
		      rcpcomm_t Comm)
{
  // read basis in a tpetra mv
  auto tp_phi = pressio::rom::test::tpetra::readBasis("basis_bdf2.txt", romSize,
						      fomSize, Comm, appObj.getDataMap());

  // construct a tpetra_block mv from a tpetra mv
  using tpb_mv = Tpetra::BlockMultiVector<>;
  tpb_mv tpb_phi(*tp_phi.data(), *appObj.getDataMap(), 1 /*1 is block size*/);
  typename decoder_d_t::jacobian_type phi(tpb_phi);

  decoder_d_t decoderObj(phi);
  return decoderObj;
}

template <typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dTpetraBlock & appObj,
		     const y1_t & yFinal, const y2_t & trueY,
		     int rank)
{
  std::string checkStr{"PASSED"};
  auto yFF_v = const_cast<y1_t &>(yFinal).data()->getVectorView().getData();
  int shift = (rank == 0) ? 0 : 10;
  const int myn = yFinal.data()->getMap()->getNodeNumElements();
  for(auto i = 0; i < myn; i++) {
    if(std::abs(yFF_v[i] - trueY[i + shift]) > 1e-10)
      checkStr = "FAILED";
    if(std::isnan(yFF_v[i]))
      checkStr = "FAILED";
  }

  return checkStr;
}

}}}//end namespace pressio::testing::wls
#endif
