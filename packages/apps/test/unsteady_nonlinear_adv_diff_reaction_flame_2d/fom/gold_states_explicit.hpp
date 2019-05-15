
#ifndef ROMPP_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_EXPLICIT_HPP_
#define ROMPP_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_EXPLICIT_HPP_

#include "ODE_ALL"

namespace rompp { namespace apps{ namespace test{

template <ode::ExplicitEnum>
struct NonLinAdvDiffReacFlame2dExpGoldStates;

template <>
struct NonLinAdvDiffReacFlame2dExpGoldStates<ode::ExplicitEnum::RungeKutta4>{
  using result_t = std::vector<double>;

  static result_t get(int Nx, int Ny,
		      double dt, double final_t){
    if(Nx==11 and Ny==21 and dt==0.001 and final_t==0.5){
      return{

      };//end return
    }else{
      std::cout << "returning empty true solution" << std::endl;
      return {};
    }
  }//get
};//struct

}}}//end namespace rompp::apps::test
#endif
