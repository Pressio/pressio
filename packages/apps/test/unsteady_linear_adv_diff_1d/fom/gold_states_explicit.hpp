
#ifndef ROMPP_DEMO_APPS_TEST_UNSTEADYLINADVDIFF1D_GOLD_EXPLICIT_HPP_
#define ROMPP_DEMO_APPS_TEST_UNTSEADYLINADVDIFF1D_GOLD_EXPLICIT_HPP_

#include "ODE_ALL"
#include "APPS_UNSTEADYLINADVDIFF1D"

namespace rompp { namespace apps{ namespace test{

template <ode::ExplicitEnum>
struct UnsteadyLinAdvDiff1dExpGoldStates;

template <>
struct UnsteadyLinAdvDiff1dExpGoldStates<ode::ExplicitEnum::Euler>{
  using result_t = std::vector<double>;
  using scalar_t = rompp::apps::UnsteadyLinAdvDiff1dEpetra::scalar_type;
  static result_t get(std::vector<scalar_t> mu, std::vector<scalar_t> domain,
		      std::vector<scalar_t> bc1D, double dt, double final_t){
    if (mu[0] == -1 && mu[1] == 1 && mu[2] == 1 && domain[0] ==0 &&
	domain[1] == 2.0 && domain[2] == 0.2 && bc1D[0] == 0 && bc1D[1] == 0 &&
	dt == 0.01 && final_t == 1.0){
      return{
	0.815361934280270,
	  1.786600784153030,
	  2.903325129361352,
	  4.123462662180391,
	  5.353273389229428,
	  6.418735933366621,
	  7.025230877624346,
	  6.701401153556792,
	  4.721667431335943};
    }else
      return{};

  } //end get
}; //end struct

}}} //end namespace rompp:apps::test
#endif
