#include "pressio/ops.hpp"
#include "tpql_operators.hpp"



class MyApp
{
public:
  MyApp(){}

  
  Eigen::Matrix<double,-1,1> createVelocity() const{
    velocity_t f(N_);
    return f;
  }
  
  template <typename state_t, typename velocity_t> 
  void velocity(const state_t & y, const double t, velocity_t & f)
  const {
    f(0) = y(0)*y(1)*y(1);
    f(1) = y(1)*y(0) + y(0)*y(0);
  }
  
  Eigen::Matrix<double,-1,1> initialCondition() const{
    state_t y(N_);
    y(0) = 3.;
    y(1) = 2.;
    return y;
  }


private: 
  int N_ = 2;
  using velocity_t = Eigen::Matrix<double,-1,1>;
  using state_t = Eigen::Matrix<double,-1,1>;

};


int main(int argc, char *argv[])
{
  using scalar_t = double;


  auto appObj = MyApp();

  int N = 2;
 
  Eigen::Matrix<scalar_t, -1,-1> Phi(N,N);
  Phi.setZero();
  Phi(0,0) = 1;
  Phi(1,1) = 1;
  auto u = appObj.initialCondition();
 
  Eigen::Matrix<scalar_t , -1,1> mu(1);
  auto PhiTJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasis(appObj, u, mu, 0., Phi);

  Eigen::Matrix<scalar_t, -1,-1> PhiTJPhi_exact(N,N);
  PhiTJPhi_exact(0,0) = 4;
  PhiTJPhi_exact(0,1) = 12;
  PhiTJPhi_exact(1,0) = 8;
  PhiTJPhi_exact(1,1) = 3;


  auto PhiTHJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasis(appObj, u, mu, 0., Phi);

  Eigen::Matrix<scalar_t, -1, 1> Hexact_list(8);
  Hexact_list(0) = 0;
  Hexact_list(1) = 4;
  Hexact_list(2) = 4;
  Hexact_list(3) = 6;
  Hexact_list(4) = 2;
  Hexact_list(5) = 1;
  Hexact_list(6) = 1;
  Hexact_list(7) = 0;

  int indx = 0;
  for (int i=0; i < 2; i++){
    for (int j=0; j < 2; j++){
      for (int k=0; k < 2 ; k++){
        std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTHJPhi[i][j][k] - Hexact_list(indx) << std::endl;
        indx += 1;}}}


  auto PhiTParamJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobian(appObj, u, mu, 0., Phi);
  std::cout << PhiTParamJPhi << std::endl;

  auto PhiTParamHPhi = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessian(appObj, u, mu, 0., Phi);
  indx = 0;
  for (int i=0; i < 2; i++){
    for (int j=0; j < 1; j++){
      for (int k=0; k < 1 ; k++){
        std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTParamHPhi[i][j][k]  << std::endl;
        indx += 1;}}}
  return 0;
}
