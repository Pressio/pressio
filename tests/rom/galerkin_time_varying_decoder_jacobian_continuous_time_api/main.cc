
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  mutable int updJacCount_ = 0;
  mutable int applyMappingCount_ = 0;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : romSize_{romSize}, jac_(fomSize, romSize)
  {
    //initialize the jacobian to be [1,2,3]
    for (int i=0; i<fomSize; ++i){
      jac_(i,0) = 1; jac_(i,1) = 2; jac_(i,2) = 3;
      jac_(i,3) = 4; jac_(i,4) = 5;
    }
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) const
  {
    ++updJacCount_;
    // updateJacobian should be called BEFORE the applyMapping.
    // Here we can use counters to verify this
    if (updJacCount_ != applyMappingCount_+1){
      throw std::runtime_error
	("invalid logic for Galerkin calls tocustomDecoder, something is off");
    }

    // at step 2, change the jacobian
    if(updJacCount_ == 2){
      for (int i=0; i<jac_.extent(0); ++i){
    	jac_(i,0) = 2; jac_(i,1) = 3; jac_(i,2) = 4;
    	jac_(i,3) = 5; jac_(i,4) = 6;
      }
    }
    // at step 3, change the jacobian again
    if(updJacCount_ == 3){
      for (int i=0; i<jac_.extent(0); ++i){
    	jac_(i,0) = 3; jac_(i,1) = 4; jac_(i,2) = 5;
    	jac_(i,3) = 6; jac_(i,4) = 7;
      }
    }
  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    ::pressio::containers::Vector<Eigen::VectorXd> & result) const
  {
    ++applyMappingCount_;
    // applyMapping should be called BEFORE the updateJacobian
    // because pressio first reconstrct FOM, then call FOM velocity
    // and then project the velocity.
    // Here we can use counters to verify this
    if (updJacCount_ != applyMappingCount_){
      throw std::runtime_error
	("invalid logic for Galerkin calls tocustomDecoder, something is off");
    }

    const auto & jacNativeObj = *jac_.data();
    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = jacNativeObj * romStateNativeObj;
  }

  const jacobian_type & getReferenceToJacobian() const{
    return jac_;
  }
};


struct MyFakeApp
{
  int N_;
public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N) : N_(N){}

  velocity_type createVelocity() const{
    return state_type(N_);
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    for (int i=0; i<N_; ++i){
      f(i) = 1;
    }
  }
};

int main(int argc, char *argv[])
{
  /* In this test we integrate: dx/dt = phi^T f( phi x)
     using Euler forward where x is vec of generalized coords.

     We then have:  x_n+1 = x_n + dt * phi^T f()

     This test is meant to test Galerkin with euler forward
     works correctly when using a decoder where its Jacobian (J_d) changes.
     To do this, we craft a test where:

     dt = 0.1, we do 3 steps so we go from t_0->t_1->t_2->t_3

     x_0 = [0 0 0 0 0]
     J_d = [1 2 3 4 5;
	    1 2 3 4 5;
	    ...
	    1 2 3 4 5];

     f=[1 1 ... 1]^T always, even if in a real case f would change.

     At step=2 (i.e. from t_1 to t_2) we update J_d to be:
     J_d = [2 3 4 5 6;
	    2 3 4 5 6;
	    ...
	    2 3 4 5 6];

     At step=3 (i.e. from t_2 to t_3) we update J_d to be:
     J_d = [3 4 5 6 7;
	    3 4 5 6 7;
	    ...
	    3 4 5 6 7]

     This means that if things are correct, we should have:
     x(t0) =  [0 0 0 0 0 ]
     x(t1) = [0 0 0 0 0 ] + 0.1 * J_d(t0)^T [1 ... 1]^T  = [1 2 3 4 5]
     x(t2) = [1 2 3 4 5 ] + 0.1 * J_d(t1)^T [1 ... 1]^T  = [3 5 7 9 11]
     x(t3) = [3 5 7 9 11] + 0.1 * J_d(t2)^T [1 ... 1]^T  = [6 9 12 15 18]
   */

  std::string checkStr {"PASSED"};

  using fom_t		= MyFakeApp;
  using scalar_t	= typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;

  using rom_state_t	= pressio::containers::Vector<Eigen::VectorXd>;
  using decoder_t	= MyCustomDecoder;

  constexpr int fomSize = 10;
  constexpr int romSize = 5;

  // app object
  fom_t appObj(fomSize);

  // decoder (use my custom one)
  decoder_t  decoderObj(fomSize, romSize);

  // this is my reference state, zero for now
  native_state_t refState(fomSize);

  // define ROM state
  rom_state_t romState(romSize);
  pressio::ops::fill(romState, 0.0);

  auto t0 = static_cast<scalar_t>(0);
  using ode_tag = pressio::ode::explicitmethods::Euler;
  using problem_t  = pressio::rom::galerkin::composeDefaultProblem<
    ode_tag, fom_t, rom_state_t, decoder_t>::type;
  problem_t galerkinProb(appObj, refState, decoderObj, romState);

  pressio::ode::advanceNSteps(galerkinProb.getStepperRef(), romState, 0.0, 0.1, 3);

  std::vector<scalar_t> trueS{6,9,12,15,18};
  for  (int i=0; i<5; ++i){
    if (std::abs(trueS[i] - romState(i)) > 1e-13) checkStr="FAILED";
  }
  std::cout << *romState.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  return 0;
}
