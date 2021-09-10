#include "pressio_rom.hpp"
#include "pressio_rom_wls.hpp"
#include "../helpers.hpp"

/*
This regression test validates the Hessian and gradient computed by the WLS system 
is correct. We test this by constructing a fake application from which we can analytically 
calculate the Hessians and gradients. 

This regression test checks BDF1 and BDF2 for upper and lower symmetric Hessian structures
*/

constexpr double dt = 1.0;

struct MyFakeApp
{
  int N_ = {};
  std::string & sentinel_;
public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N, std::string & sentinel)
    : N_(N), sentinel_(sentinel){}

  velocity_type createVelocity() const{
    return state_type(N_);
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const
  {
    return dense_matrix_type(N_, B.cols());
  }

  void applyJacobian(const state_type &,
  		     const dense_matrix_type & B,
  		     scalar_type time,
  		     dense_matrix_type & A) const
  {
    A.setConstant(time + 1.);
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    f.setConstant(time + 1.);
  }
};

template<typename gradient_t, typename hessian_t>
struct MyFakeSolver
{
  int fomSize_ = {};
  int romSize_ = {};
  int callCounter_ = 0;
  gradient_t gradient_;
  hessian_t H_;
  std::string & checkString_;
  Eigen::MatrixXd trueH1_;
  Eigen::MatrixXd trueH2_;
  Eigen::VectorXd trueGradient1_;
  Eigen::VectorXd trueGradient2_;
  //Constructor for BDF1 with lower triangular Hessian
  MyFakeSolver(int fomSize, int romSize,
	       std::string & checkString, 
               ::pressio::ode::BDF1 odeTag,
               ::pressio::matrixLowerTriangular hessianStructureTag)
    : fomSize_(fomSize),
      romSize_(romSize),
      gradient_(romSize),
      H_(romSize, romSize),
      checkString_(checkString)
  {
    trueH1_.resize(romSize,romSize);
    trueH2_.resize(romSize,romSize);
    trueH1_.setZero(romSize,romSize);
    trueH2_.setZero(romSize,romSize);
    trueGradient1_.resize(romSize);
    trueGradient2_.resize(romSize);
    trueGradient1_.setZero();
    trueGradient2_.setZero();
      //only fill lower block
      for (auto i=0; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(i,j) = 3.;}}
      for (auto i=2; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(i,j) = 6.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH1_(i,j) = 12.;}}
      trueGradient1_(0) = 6.; trueGradient1_(1) = 6.; trueGradient1_(2) = 15.; 
      trueGradient1_(3) = 15.; trueGradient1_(4) = 18.; trueGradient1_(5) = 18.;

      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(i,j) = 30.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(i,j) = 12.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(i,j) = 51.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(i,j) = 15.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH2_(i,j) = 75.;}}

      trueGradient2_(0) = 51.; trueGradient2_(1) = 51.; trueGradient2_(2) = 78.; 
      trueGradient2_(3) = 78.; trueGradient2_(4) = 90.; trueGradient2_(5) = 90.;

  }

  //Constructor for BDF1 with upper triangular Hessian
  MyFakeSolver(int fomSize, int romSize,
	       std::string & checkString, 
               ::pressio::ode::BDF1 odeTag,
               ::pressio::matrixUpperTriangular hessianStructureTag)
    : fomSize_(fomSize),
      romSize_(romSize),
      gradient_(romSize),
      H_(romSize, romSize),
      checkString_(checkString)
  {
    trueH1_.resize(romSize,romSize);
    trueH2_.resize(romSize,romSize);
    trueH1_.setZero(romSize,romSize);
    trueH2_.setZero(romSize,romSize);
    trueGradient1_.resize(romSize);
    trueGradient2_.resize(romSize);
    trueGradient1_.setZero();
    trueGradient2_.setZero();
      //only fill upper block
      for (auto i=0; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(j,i) = 3.;}}
      for (auto i=2; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(j,i) = 6.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH1_(j,i) = 12.;}}
      trueGradient1_(0) = 6.; trueGradient1_(1) = 6.; trueGradient1_(2) = 15.; 
      trueGradient1_(3) = 15.; trueGradient1_(4) = 18.; trueGradient1_(5) = 18.;

      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(j,i) = 30.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(j,i) = 12.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(j,i) = 51.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(j,i) = 15.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH2_(j,i) = 75.;}}

      trueGradient2_(0) = 51.; trueGradient2_(1) = 51.; trueGradient2_(2) = 78.; 
      trueGradient2_(3) = 78.; trueGradient2_(4) = 90.; trueGradient2_(5) = 90.;

  }


  //Constructor for BDF2 with lower triangular Hessian
  MyFakeSolver(int fomSize, int romSize,
	       std::string & checkString, 
               ::pressio::ode::ode::BDF2 odeTag,
               ::pressio::matrixLowerTriangular hessianStructureTag)
    : fomSize_(fomSize),
      romSize_(romSize),
      gradient_(romSize),
      H_(romSize, romSize),
      checkString_(checkString)
  {
    trueH1_.resize(romSize,romSize);
    trueH2_.resize(romSize,romSize);
    trueH1_.setZero(romSize,romSize);
    trueH2_.setZero(romSize,romSize);
    trueGradient1_.resize(romSize);
    trueGradient2_.resize(romSize);
    trueGradient1_.setZero();
    trueGradient2_.setZero();
      //only fill lower block
      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(i,j) = 5. + 2./3.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(i,j) = 5. + 2./3.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(i,j) = -1.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(i,j) = 4.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH1_(i,j) = 3.;}}

      trueGradient1_(0) = 5.*2./3.; trueGradient1_(1) = 5.*2./3.; trueGradient1_(2) = 14.*2./3.; 
      trueGradient1_(3) = 14.*2./3.; trueGradient1_(4) = 9.*2./3.; trueGradient1_(5) = 9.*2./3.;

      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(i,j) = 14.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(i,j) = 8.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(i,j) = -3.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(i,j) = 21. + 2./3.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(i,j) = 12.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH2_(i,j) = 27.;}}

      trueGradient2_(0) = 34.*2./3.; trueGradient2_(1) = 34.*2./3.; trueGradient2_(2) = 59.*2./3.; 
      trueGradient2_(3) = 59.*2./3.; trueGradient2_(4) = 54.*2./3.; trueGradient2_(5) = 54.*2./3.;
  }

  //Constructor for BDF2 with upper triangular Hessian
  MyFakeSolver(int fomSize, int romSize,
	       std::string & checkString, 
               ::pressio::ode::ode::BDF2 odeTag,
               ::pressio::matrixUpperTriangular hessianStructureTag)
    : fomSize_(fomSize),
      romSize_(romSize),
      gradient_(romSize),
      H_(romSize, romSize),
      checkString_(checkString)
  {
    trueH1_.resize(romSize,romSize);
    trueH2_.resize(romSize,romSize);
    trueH1_.setZero(romSize,romSize);
    trueH2_.setZero(romSize,romSize);
    trueGradient1_.resize(romSize);
    trueGradient2_.resize(romSize);
    trueGradient1_.setZero();
    trueGradient2_.setZero();

      //only fill upper block
      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(j,i) = 5. + 2./3.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(j,i) = 5. + 2./3.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=0; j<2; ++j){
          trueH1_(j,i) = -1.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH1_(j,i) = 4.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH1_(j,i) = 3.;}}

      trueGradient1_(0) = 5.*2./3.; trueGradient1_(1) = 5.*2./3.; trueGradient1_(2) = 14.*2./3.; 
      trueGradient1_(3) = 14.*2./3.; trueGradient1_(4) = 9.*2./3.; trueGradient1_(5) = 9.*2./3.;

      for (auto i=0; i<2; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(j,i) = 14.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(j,i) = 8.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=0; j<2; ++j){
          trueH2_(j,i) = -3.;}}
      for (auto i=2; i<4; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(j,i) = 21. + 2./3.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=2; j<4; ++j){
          trueH2_(j,i) = 12.;}}
      for (auto i=4; i<6; ++i){
        for (auto j=4; j<6; ++j){
          trueH2_(j,i) = 27.;}}

      trueGradient2_(0) = 34.*2./3.; trueGradient2_(1) = 34.*2./3.; trueGradient2_(2) = 59.*2./3.; 
      trueGradient2_(3) = 59.*2./3.; trueGradient2_(4) = 54.*2./3.; trueGradient2_(5) = 54.*2./3.;
  }

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++callCounter_;
    double rnorm = 0.;
    bool recomputeJacobian = true;
    sys.hessianAndGradient(state,H_,gradient_,::pressio::Norm::L2,rnorm,recomputeJacobian);
    auto native_state = *state.data();
    if (callCounter_==1)
    {
      if (!trueGradient1_.isApprox(*gradient_.data())) checkString_ = "FAILED";
      if (!trueH1_.isApprox(*H_.data())) checkString_ = "FAILED";
    }
    if (callCounter_==2)
    {
      if (!trueGradient2_.isApprox(*gradient_.data())) checkString_ = "FAILED";
      if (!trueH2_.isApprox(*H_.data())) checkString_ = "FAILED";
    }

  }
};

template< typename ode_tag,typename hessian_matrix_structure_tag>
std::string doRun()
{
  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using native_state_t  = typename fom_t::state_type;
  using fom_state_t     = pressio::containers::Vector<native_state_t>;
  using scalar_t        = typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using wls_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_t	= MyCustomDecoder;
  using wls_hessian_t	= typename ::pressio::containers::DenseMatrix<Eigen::Matrix<scalar_t, -1, -1>>;

  constexpr int fomSize = 3;
  constexpr int romSize = 2;
  constexpr int numStepsInWindow = 3;
  constexpr int numWindows = 2;

  // app object
  fom_t appObj(fomSize, checkStr);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);
  refState.setConstant(0.0);
  fom_state_t refStateWrapper(refState);

  // define ROM state
  wls_state_t wlsState(romSize*numStepsInWindow);
  wls_state_t wlsStateIC(romSize);
  pressio::ops::fill(wlsState, 0.0);
  pressio::ops::fill(wlsStateIC, 1.0);

  //*** WLS problem ***
  using precon_type = ::pressio::rom::wls::preconditioners::NoPreconditioner;
  using jacobians_update_tag = ::pressio::rom::wls::NonFrozenJacobian;
  using policy_t     = pressio::rom::wls::HessianGradientSequentialPolicy
    <fom_t, decoder_t, ode_tag, hessian_matrix_structure_tag, precon_type, jacobians_update_tag>;
  using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi
    <wls_state_t, decoder_t, wls_hessian_t, policy_t>;

  policy_t hgPolicy(romSize, numStepsInWindow, decoderObj,
		    appObj, refStateWrapper);

 
  wls_system_t wlsSystem(decoderObj, hgPolicy, refStateWrapper,
                         refStateWrapper, wlsStateIC);


  using solver_t = MyFakeSolver<wls_state_t,wls_hessian_t>;
  solver_t solver(fomSize*numStepsInWindow, romSize*numStepsInWindow, checkStr,ode_tag(),hessian_matrix_structure_tag());

  ::pressio::rom::wls::solveWindowsSequentially(wlsSystem, wlsState, solver, numWindows, dt);
  return checkStr;

}
int main(int argc, char *argv[])
{
  std::string checkStr = "PASSED";

  // Test for BDF1 with lower triangular
  auto checkStr1 = doRun<::pressio::ode::BDF1,
                          ::pressio::matrixLowerTriangular>();
  if (checkStr1 == "FAILED"){
    std::cout << "WLS failed on implicit Euler with lower triangular Hessian" << std::endl;
    checkStr = "FAILED";
  }
  // Test for BDF1 with upper triangular
  auto checkStr2 = doRun<::pressio::ode::BDF1,
                          ::pressio::matrixUpperTriangular>();
  if (checkStr2 == "FAILED"){
    std::cout << "WLS failed on implicit Euler with upper triangular Hessian" << std::endl;
    checkStr = "FAILED";
  }

   // Test for BDF2 with lower triangular
  auto checkStr3 = doRun<::pressio::ode::ode::BDF2,
                          ::pressio::matrixLowerTriangular>();
  if (checkStr3 == "FAILED"){
     std::cout << "WLS failed on BDF2 with lower triangular Hessian" << std::endl;
     checkStr = "FAILED";
  }

   // Test for BDF2 with upper triangular
  auto checkStr4 = doRun<::pressio::ode::ode::BDF2,
                          ::pressio::matrixUpperTriangular>();
  if (checkStr4 == "FAILED"){
    std::cout << "WLS failed on BDF2 with upper triangular Hessian" << std::endl;
    checkStr = "FAILED";
  }
  std::cout << checkStr << std::endl;
  return 0;
}

