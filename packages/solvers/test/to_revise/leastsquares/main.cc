
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

// #include "algebra_ConfigDefs.hpp"
// #include "algebra_meta.hpp"
// #include "algebra_forward_declarations.hpp"
// #include "algebra_static_assert.hpp"

#include "vector/algebra_vector_eigen.hpp"
#include "matrix/algebra_matrix_traits.hpp"
#include "matrix/algebra_matrix_eigen.hpp"
#include <Eigen/Dense>
#include "least_squares/algebra_leastsquares.hpp"


double g(double x, const Eigen::VectorXd & a){
  return a[0]*x/(a[1]+x);
  //  return a[0]*std::exp(a[1]*x)*std::cos(a[2]*x+a[3]);
}
double drdth0(double x, const Eigen::VectorXd & a){
  return -x/(a[1]+x);
  //  return -std::exp(a[1]*x)*std::cos(a[2]*x+a[3]);
}
double drdth1(double x, const Eigen::VectorXd & a){
  return a[0]*x/( (a[1]+x)*(a[1]+x) );
  //  return -x*a[0]*std::exp(a[1]*x)*std::cos(a[2]*x+a[3]);
}
// double drdth2(double x, const Eigen::VectorXd & a){
//   return a[0]*std::exp(a[1]*x)*std::sin(a[2]*x+a[3])*x;
// }
// double drdth3(double x, const Eigen::VectorXd & a){
//   return a[0]*std::exp(a[1]*x)*std::sin(a[2]*x+a[3]);
// }


struct app
{
  using state_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;

  int N_ = 7;
  state_type measX;
  state_type measY;

  void init(){
    measX << 0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740;
    measY << 0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317;
    // measX = {0., 1.11111111, 2.22222222, 3.33333333,
    // 	     4.44444444, 5.55555556, 6.66666667,
    // 	     7.77777778,  8.88888889,  10.};
    // measY = {0.59866925, 0.23428037, -0.34075947, 0.22351218,
    // 	     0.01224858, -0.02007828, -0.0214535,
    // 	     -0.02288888,  0.02298912,  0.02216982};  
  }
  
  void operator() ( const state_type & u,
  		    state_type & R,
  		    jacobian_type & jac)
  {
    // R.resize(N_);
    // jac.resize(N_,2);
    // auto & jacEmat = jac.getNonConstRefToDataImpl();
    // double sum = 0.0;
    // for (int i=0; i<N_; i++)
    // {
    //  R[i] = measY[i] - g(measX[i], *u.view() );
    //  sum += R[i]*R[i];
    //  jacEmat(i,0) = drdth0(measX[i], *u.view());
    //  jacEmat(i,1) = drdth1(measX[i], *u.view());
    //  // jacEmat(i,2) = drdth2(measX[i], *u.view());
    //  // jacEmat(i,3) = drdth3(measX[i], *u.view());
    // }    
    // std::cout << "sum " << sum << std::endl;

  };//end op()

};

  
int main()
{     
  // using myvec_t = algebra::Vector<app::state_type>;
  // using mymat_t = algebra::Matrix<app::jacobian_type>;

  // app::myvec_t sol;
  // sol.resize(2);
  // sol[0] = 0.01;
  // sol[1] = 0.1;
  
  // app F;
  // F.init();
  // algebra::nonLinearLstsq<app,app::myvec_t,app::mymat_t>(F, sol);
  // std::cout << sol[0] << " "
  // 	    << sol[1]
  // 	    << std::endl;

  return 0;
};



  // Eigen::MatrixXd A;
  // A.resize(3,2);
  // A(0,0) =  0.68;
  // A(0,1) = 0.597;
  // A(1,0) = -0.211;
  // A(1,1) = 0.823;
  // A(2,0) = 0.566;
  // A(2,1) = -0.605;      
  // app::mymat_t myA(A);
  // std::cout << "Here is the matrix A:\n" << *myA.view() << std::endl;

  // Eigen::VectorXd b; b.resize(3);
  // b(0) = -0.33;
  // b(1) = 0.536;
  // b(2) = -0.444;

  // app::myvec_t myb(b);
  // std::cout << "Here is the vec b:\n" << *myb.view() << std::endl;
  // app::myvec_t myx;
  // myx.resize(myA.cols());
  // algebra::linearLstsq(myA, myb, myx);
  // std::cout << myx[0] << " " << myx[1] << std::endl;
  
