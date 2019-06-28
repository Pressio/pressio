
#include <gtest/gtest.h>
#include "ALGEBRA_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"

TEST(adv_diff_reaction_2d_eigen, spatial_residual){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;

  static_assert(std::is_same<scalar_t, double>::value, "");

  constexpr auto zero = ::rompp::algebra::constants::zero<scalar_t>();

  constexpr int Nx = 5, Ny = 5;
  app_t appobj(Nx, Ny);
  appobj.setup();

  const auto totUnk = appobj.getUnknownCount();
  EXPECT_EQ( totUnk, 45 );

  // set fake init state to test:
  // it is important to set state to be something
  // different at every dof otherwise we might be tricked
  app_state_t yTest(totUnk);
  ::rompp::algebra::Vector<app_state_t> y(yTest);
  for (auto i=0; i<y.size(); ++i)
    y[i] = i*0.1;

  //compute residual
  const auto r0n = appobj.residual(*y.data(), zero);

  const auto eps = appobj.getDiffusivity();
  const auto K = appobj.getReactionRate();
  const auto dx = appobj.getDx();
  const auto dy = appobj.getDy();
  EXPECT_DOUBLE_EQ( dx, 0.25);
  EXPECT_DOUBLE_EQ( dy, 0.5);

  // remember that u,v,s1,s2,s3
  // have values at grid points, while state and residual
  // have more values because of the dofs
  const auto u = appobj.getU();
  const auto v = appobj.getV();
  const auto s1 = appobj.getS1();
  const auto s2 = appobj.getS2();
  const auto s3 = appobj.getS3();

  constexpr scalar_t tol = 1e-14;

  int iGPt = 0;// enumerates grid points
  for (auto i=0; i<18; ++i){
    const auto A = eps/(dx*dx) + eps/(dy*dy);
    const auto B = -u(iGPt)/(2.*dx) + eps/(dx*dx);
    const auto F = +u(iGPt)/(2.*dx) + eps/(dx*dx);
    const auto D = -v(iGPt)/(2.*dy) + eps/(dy*dy);
    const auto E = +v(iGPt)/(2.*dy) + eps/(dy*dy);

    const auto appResValue = r0n(i);
    scalar_t trueValue{};

    // i=0,1,2 cover the dofs at grid point 0
    if (i==0){
      trueValue = -2.*A*y(0) + B*y(3) + (D+E)*y(9)-K*y(0)*y(1) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    if (i==1){
      trueValue = -2.*A*y(1) + B*y(4) + (D+E)*y(10)-K*y(0)*y(1) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    if (i==2){
      trueValue = -2.*A*y(2) + B*y(5) + (D+E)*y(11)+K*y(0)*y(1)
	-K*y(2) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=3,4,5 cover the dofs at grid point 1
    if (i==3){
      trueValue = -2.*A*y(3) + B*y(6) + F*y(0)
    	+(D+E)*y(12)-K*y(3)*y(4) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==4){
      trueValue = -2.*A*y(4) + B*y(7) + F*y(1)
    	+(D+E)*y(13)-K*y(3)*y(4) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==5){
      trueValue = -2.*A*y(5) + B*y(8) + F*y(2)
    	+(D+E)*y(14)+K*y(3)*y(4) -K*y(5) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=6,7,8 cover the dofs at grid point 2
    if (i==6){
      trueValue = -2.*A*y(6) + F*y(3)
    	+(D+E)*y(15)-K*y(6)*y(7) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==7){
      trueValue = -2.*A*y(7) + F*y(4)
    	+(D+E)*y(16)-K*y(6)*y(7) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==8){
      trueValue = -2.*A*y(8) + F*y(5)
    	+(D+E)*y(17)+K*y(6)*y(7) -K*y(8) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=9,10,11 cover the dofs at grid point 3
    if (i==9){
      trueValue = -2.*A*y(9) + B*y(12) +
    	+D*y(18)+E*y(0) -K*y(9)*y(10) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==10){
      trueValue = -2.*A*y(10) + B*y(13)
    	+D*y(19)+E*y(1)-K*y(9)*y(10) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==11){
      trueValue = -2.*A*y(11) + B*y(14)
    	+D*y(20)+E*y(2)+K*y(9)*y(10) -K*y(11) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=12,13,14 cover the dofs at grid point 4
    if (i==12){
      trueValue = -2.*A*y(12) + B*y(15) + F*y(9)
    	+D*y(21)+E*y(3) -K*y(12)*y(13) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==13){
      trueValue = -2.*A*y(13) + B*y(16) + F*y(10)
    	+D*y(22)+E*y(4) -K*y(12)*y(13) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==14){
      trueValue = -2.*A*y(14) + B*y(17) + F*y(11)
    	+D*y(23)+E*y(5) +K*y(12)*y(13) - K*y(14) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=36,37,38
    if (i==36){
      trueValue = -2.*A*y(36) + B*y(39)
    	+(D+E)*y(27) -K*y(36)*y(37) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==37){
      trueValue = -2.*A*y(37) + B*y(16)
    	+(D+E)*y(28) -K*y(36)*y(37) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==38){
      trueValue = -2.*A*y(38) + B*y(17)
    	+(D+E)*y(29) +K*y(36)*y(37) - K*y(38) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    // i=39,40,41
    if (i==39){
      trueValue = -2.*A*y(39) + B*y(42) + F*y(36)
    	+(D+E)*y(30) -K*y(39)*y(40) + s1(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==40){
      trueValue = -2.*A*y(40) + B*y(43) + F*y(37)
    	+(D+E)*y(31) -K*y(39)*y(40) + s2(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }
    if (i==41){
      trueValue = -2.*A*y(41) + B*y(44) + F*y(38)
    	+(D+E)*y(32) +K*y(39)*y(40) - K*y(41) + s3(iGPt);
      EXPECT_NEAR( trueValue, appResValue, tol);
    }

    if ( (i+1) % 3 ==0 )
      iGPt++;
  }
}
