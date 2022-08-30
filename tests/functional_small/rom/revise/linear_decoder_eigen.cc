
#include <gtest/gtest.h>
#include "pressio/rom_decoder.hpp"

TEST(rom, linear_decoder_eigen)
{
  using fom_state_t	  = Eigen::VectorXd;
  using decoder_jac_t	= Eigen::MatrixXd;

  decoder_jac_t A(10, 3); 
  A.fill(1.);

  auto decoder = pressio::rom::create_time_invariant_linear_decoder<fom_state_t>(A);

  // construct reduced state
  using rom_state_t = Eigen::VectorXd;
  rom_state_t yRom(A.cols());
  yRom.setConstant(2.);

  // apply mapping
  fom_state_t yFom(A.rows());
  decoder.applyMapping(yRom, yFom);

  for (auto i=0; i<yFom.size(); ++i){
    EXPECT_DOUBLE_EQ(yFom(i), 6.);
  }
}
