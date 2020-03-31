#ifndef BURGERS_FOM_FUNCTIONS_KOKKOS_ 

#include "utils_kokkos.hpp"
namespace{
//for Kokkos=================================
template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dKokkos & appObj, ::pressio::ode::implicitmethods::Euler & odeTag, int  romSize, int  fomSize){

  using fom_dmat_t = typename decoder_d_t::jacobian_type;
  fom_dmat_t phi("phi", fomSize, 11);
  pressio::rom::test::kokkos::readBasis("basis_euler.txt", romSize, fomSize, *phi.data());
  decoder_d_t decoderObj(phi);

  return decoderObj;
}

template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dKokkos & appObj, ::pressio::ode::implicitmethods::BDF2 & odeTag, int  romSize, int  fomSize){

  using fom_dmat_t = typename decoder_d_t::jacobian_type;
  fom_dmat_t phi("phi", fomSize, 11);
  pressio::rom::test::kokkos::readBasis("basis_bdf2.txt", romSize, fomSize, *phi.data());
  decoder_d_t decoderObj(phi);

  return decoderObj;
}

template<typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dKokkos & appObj ,y1_t yFinal,y2_t trueY, int fomSize){

  std::string checkStr {"PASSED"};

  using native_state_t_h = pressio::apps::Burgers1dKokkos::state_type_h;
  native_state_t_h yFinal_h("yFF_h", fomSize);
  Kokkos::deep_copy(yFinal_h, *yFinal.data());

  for (auto i=0; i<fomSize; i++){
    if (std::abs(yFinal_h(i) - trueY[i]) > 1e-8) checkStr = "FAILED";
  }

  return checkStr;
}
//==============================================
}//end namespace

#endif
