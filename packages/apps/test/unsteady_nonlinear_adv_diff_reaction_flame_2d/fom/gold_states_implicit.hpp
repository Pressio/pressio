
#ifndef ROMPP_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_IMPLICIT_HPP_
#define ROMPP_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_IMPLICIT_HPP_

#include "ODE_ALL"

namespace rompp { namespace apps{ namespace test{

template <ode::ImplicitEnum>
struct NonLinAdvDiffReacFlame2dImpGoldStates;

template <>
struct NonLinAdvDiffReacFlame2dImpGoldStates<ode::ImplicitEnum::Euler>{
  using result_t = std::vector<double>;

  static result_t get(int Nx, int Ny,
		      double dt, double final_t){
    if(Nx==11 and Ny==21 and dt==0.1 and final_t==1.0){
      return{

      };//end return
    }else{
      std::cout << "returning empty true solution" << std::endl;
      return {};
    }
  }//get
};//struct


template <>
struct NonLinAdvDiffReacFlame2dImpGoldStates<ode::ImplicitEnum::BDF2>{
  using result_t = std::vector<double>;

  static result_t get(int Nx, int Ny,
		      double dt, double final_t){
    if(Nx==11 and Ny==21 and dt==0.1 and final_t==1.0){
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
