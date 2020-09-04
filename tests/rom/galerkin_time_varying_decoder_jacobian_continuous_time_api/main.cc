
#include "pressio_rom.hpp"
#include "pressio_apps.hpp"
#include "utils_eigen.hpp"
#include "custom_mapping.hpp"

struct MyFakeApp
{
  int N_;
public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;

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
  problem_t galerkinProb(appObj, refState, decoderObj, romState, t0);

  pressio::ode::advanceNSteps(galerkinProb.getStepperRef(), romState, 0.0, 0.1, 3);

  std::vector<scalar_t> trueS{6,9,12,15,18};
  for  (int i=0; i<5; ++i){
    if (std::abs(trueS[i] - romState(i)) > 1e-13) checkStr="FAILED";
  }
  std::cout << *romState.data() << std::endl;

  std::cout << checkStr <<  std::endl;
  return 0;
}
