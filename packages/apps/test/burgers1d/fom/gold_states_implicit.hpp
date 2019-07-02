
#ifndef PRESSIO_DEMO_APPS_TEST_BURGERS1D_GOLD_IMPLICIT_HPP_
#define PRESSIO_DEMO_APPS_TEST_BURGERS1D_GOLD_IMPLICIT_HPP_

#include "ODE_ALL"

namespace pressio { namespace apps{ namespace test{

template <ode::ImplicitEnum>
struct Burgers1dImpGoldStates;

template <>
struct Burgers1dImpGoldStates<ode::ImplicitEnum::Euler>{
  using result_t = std::vector<double>;

  static result_t get(int N,
		      double dt,
		      double final_t = 35.,
		      double mu1 = 5.0,
		      double mu2 = 0.02,
		      double mu3 = 0.02)
  {
    if (N == 20 and dt == 0.01 and final_t == 35. and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return { 5.0209814000128, 5.044067908724,
	       5.0694601439534, 5.0973757621517,
	       5.1280505160622, 5.1617393078782,
	       5.1987172219743, 5.2392805124856,
	       5.2837474958645, 5.33245922794,
	       5.3857796605719, 5.4440945290029,
	       5.507807234196,  5.577327956428,
	       5.6530483940524, 5.7352878834556,
	       5.8241864748671, 5.919507584836,
	       6.0203022968451, 6.1243942597171};
    }
    else if (N == 20 and dt == 0.01 and final_t == 0.10 and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return {1.2392405345107, 1.0051378268469,
	      1.0025875046782, 1.0028353031206,
	      1.0031333374311, 1.0034628717396,
	      1.0038270641633, 1.0042295588277,
	      1.0046743839626, 1.0051659914443,
	      1.0057093013441, 1.0063097511659,
	      1.0069733502617, 1.0077067399692,
	      1.0085172600729, 1.0094130222541,
	      1.0104029912645, 1.0114970746344,
	      1.0127062218147, 1.0140425337419};
    }
    else
      return {};

  }//end get
};//end struct


template <>
struct Burgers1dImpGoldStates<ode::ImplicitEnum::BDF2>{
  using result_t = std::vector<double>;

  static result_t get(int N,
		      double dt,
		      double final_t = 35.,
		      double mu1 = 5.0,
		      double mu2 = 0.02,
		      double mu3 = 0.02)
  {
    if (N == 20 and dt == 0.01 and final_t == 20. and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return {5.0206895255691, 5.0404333479369,
	      5.0457472457546, 4.9888777834243,
	      4.7415195052177, 4.0783524696869,
	      2.9570513683612, 1.9628122693434,
	      1.5544101171779, 1.4851802719293,
	      1.5102650056944, 1.5578826205313,
	      1.6138128119771, 1.6758489014741,
	      1.7440245069426, 1.8187984884494,
	      1.9007347577867, 1.9904481702812,
	      2.0885961580656, 2.1958778186853};
    }
    else if (N == 20 and dt == 0.01 and final_t == 0.10 and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return {1.2394672849251, 1.004917457916,
	      1.0025833517661, 1.0028354862953,
	      1.0031335910993, 1.0034631525729,
	      1.0038273746478, 1.0042299021045,
	      1.0046747635112, 1.0051664111167,
	      1.005709765406, 1.0063102643415,
	      1.0069739177846, 1.0077073676384,
	      1.0085179543157, 1.0094137901963,
	      1.0104038408089, 1.0114980145485,
	      1.0127072618291, 1.0140436846607};
    }
    else
      return {};

  }//end get
};//end struct

}}}//end namespace pressio::apps::test
#endif
