
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"

TEST(adv_diff_reaction_2d_eigen, spatial_jacobian){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;

  static_assert(std::is_same<scalar_t, double>::value, "");

  constexpr auto one = ::rompp::core::constants::one<scalar_t>();
  constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

  constexpr int Nx = 5, Ny = 5;
  app_t appobj(Nx, Ny);
  appobj.setup();

  const auto totUnk = appobj.getUnknownCount();
  EXPECT_EQ( totUnk, 45 );

  // set fake init state to test:
  // it is important to set state to be something
  // different at every dof otherwise we might be tricked
  app_state_t yTest(totUnk);
  ::rompp::core::Vector<app_state_t> y(yTest);
  for (auto i=0; i<y.size(); ++i)
    y[i] = i*0.1;

  //compute jacobian (sparse)
  const auto Js = appobj.jacobian(*y.data(), zero);
  // convert from sparse to dense Eigen matrix
  // this is for convenience for testing purposes
  const Eigen::MatrixXd Jd(Js);

  std::cout << Js << std::endl;

  std::cout << Js.coeff(0,0) << std::endl;

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
  for (auto i=0; i<45; ++i){
    const auto A = eps/(dx*dx) + eps/(dy*dy);
    const auto B = -u(iGPt)/(2.*dx) + eps/(dx*dx);
    const auto F = +u(iGPt)/(2.*dx) + eps/(dx*dx);
    const auto D = -v(iGPt)/(2.*dy) + eps/(dy*dy);
    const auto E = +v(iGPt)/(2.*dy) + eps/(dy*dy);

    scalar_t trueValue{};

    // i=0,1,2 cover the dofs at grid point 0
    if (i==0){
      EXPECT_NEAR(-2.*A-K*y(1), Js.coeff(i,i), tol);
      EXPECT_NEAR( B, Js.coeff(i,3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,9), tol);
      EXPECT_NEAR( -K*y(0), Js.coeff(i,1), tol);
    }
    if (i==1){
      EXPECT_NEAR(-2.*A-K*y(0), Js.coeff(i,i), tol);
      EXPECT_NEAR( B, Js.coeff(i,4), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,10), tol);
      EXPECT_NEAR( -K*y(1), Js.coeff(i,0), tol);
    }
    if (i==2){
      EXPECT_NEAR(-2.*A-K, Js.coeff(i,i), tol);
      EXPECT_NEAR( B, Js.coeff(i,5), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,11), tol);
      EXPECT_NEAR( K*y(i-1), Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2), Js.coeff(i,i-1), tol);
    }

    // i=3,4,5 cover the dofs at grid point 1
    if (i==3){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==4){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==5){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( B,		Js.coeff(i,i+3), tol);
      EXPECT_NEAR( F,		Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E,		Js.coeff(i,i+9), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=6,7,8 cover the dofs at grid point 2
    if (i==6){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==7){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==8){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( F,		Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D+E,		Js.coeff(i,i+9), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=9,10,11 cover the dofs at grid point 3
    if (i==9 or i==10 or i==11){
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( D, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( E, Js.coeff(i,i-9), tol);
    }
    if (i==9){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==10){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==11){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=12,13,14 cover the dofs at grid point 4
    if (i==12 or i==13 or i==14){
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( E, Js.coeff(i,i-9), tol);
    }
    if (i==12){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==13){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==14){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=15,16,17 cover the dofs at grid point 5
    if (i==15 or i==16 or i==17){
      EXPECT_NEAR( F, Js.coeff(i,i-3), tol);
      EXPECT_NEAR( D, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( E, Js.coeff(i,i-9), tol);
    }
    if (i==15){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==16){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==17){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=27,28,29 cover the dofs at grid point 9
    if (i==27 or i==28 or i==29){
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( D, Js.coeff(i,i+9), tol);
      EXPECT_NEAR( E, Js.coeff(i,i-9), tol);
    }
    if (i==27){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==28){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==29){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    // i=36,37,38 cover the dofs at grid point 12
    if (i==36 or i==37 or i==38){
      EXPECT_NEAR( B, Js.coeff(i,i+3), tol);
      EXPECT_NEAR( D+E, Js.coeff(i,i-9), tol);
    }
    if (i==36){
      EXPECT_NEAR(-2.*A-K*y(i+1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i+1), tol);
    }
    if (i==37){
      EXPECT_NEAR(-2.*A-K*y(i-1), Js.coeff(i,i), tol);
      EXPECT_NEAR( -K*y(i), Js.coeff(i,i-1), tol);
    }
    if (i==38){
      EXPECT_NEAR(-2.*A-K,	Js.coeff(i,i), tol);
      EXPECT_NEAR( K*y(i-1),	Js.coeff(i,i-2), tol);
      EXPECT_NEAR( K*y(i-2),	Js.coeff(i,i-1), tol);
    }

    if ( (i+1) % 3 ==0 )
      iGPt++;
  }
}
