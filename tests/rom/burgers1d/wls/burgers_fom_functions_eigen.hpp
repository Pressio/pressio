#include "utils_eigen.hpp"
namespace{

template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dEigen & appObj, ::pressio::ode::implicitmethods::Euler & odeTag, int  romSize, int  fomSize){

  const auto phiNative = pressio::rom::test::eigen::readBasis("basis_euler.txt", romSize, fomSize);
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

template <typename decoder_d_t>
decoder_d_t readBasis( pressio::apps::Burgers1dEigen & appObj, ::pressio::ode::implicitmethods::BDF2 & odeTag, int  romSize, int  fomSize){

  const auto phiNative = pressio::rom::test::eigen::readBasis("basis_bdf2.txt", romSize, fomSize);
  decoder_d_t decoderObj(phiNative);

  return decoderObj;
}

template<typename y1_t, typename y2_t>
std::string checkSol(pressio::apps::Burgers1dEigen & appObj ,y1_t yFinal,y2_t trueY, int fomSize){

  std::string checkStr {"PASSED"};

  for (int i=0;i<fomSize;i++){
    if (std::abs(yFinal[i] - trueY[i]) > 1e-10) checkStr = "FAILED";
  }
  return checkStr;
}

} //end namespace
