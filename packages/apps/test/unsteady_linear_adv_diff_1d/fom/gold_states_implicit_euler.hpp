
#ifndef ROMPP_APPS_TEST_UNSTEADY_LIN_ADV_DIFF1D_GOLD_IMPLICIT_EULER_HPP_
#define ROMPP_APPS_TEST_UNSTEADY_LIN_ADV_DIFF1D_GOLD_IMPLICIT_EULER_HPP_

namespace rompp { namespace apps{ namespace test{

template <ode::ImplicitEnum>
struct UnsteadyLinAdvDiff1dImplicitGoldStates;

template<>
struct UnsteadyLinAdvDiff1dImplicitGoldStates<ode::ImplicitEnum::Euler>{
  using result_t = std::vector<double>;
  using scalar_t = rompp::apps::UnsteadyLinAdvDiff1dEpetra::scalar_type;
  static result_t get(std::vector<scalar_t> mu, std::vector<scalar_t>
		      domain, std::vector<scalar_t> bc1D, double dt,
		      double final_t){
    if (mu[0] == -0.857241631161166 &&
	mu[1] ==  0.104833925269630 &&
	mu[2] == -0.713183149274631 &&
	domain[0] ==0 && domain[1] == 2.0 && domain[2] == 0.05 &&
	bc1D[0] == 0 && bc1D[1] == 0.25 &&
        dt == 0.1 && final_t == 5.0){
      return{
	     0.078785663735240,
	     0.157624101933472,
	     0.236101732342969,
	     0.313818822774520,
	     0.390389396123673,
	     0.465441111166957,
	     0.538615121426790,
	     0.609565914251985,
	     0.677961132121064,
	     0.743481378043633,
	     0.805820006810428,
	     0.864682903724941,
	     0.919788252338458,
	     0.970866292605498,
	     1.017659070777805,
	     1.059920182261742,
	     1.097414508576110,
	     1.129917949464554,
	     1.157217151138739,
	     1.179109231555018,
	     1.195401503558213,
	     1.205911196661081,
	     1.210465178166923,
	     1.208899674285308,
	     1.201059991836914,
	     1.186800241092804,
	     1.165983060245900,
	     1.138479341967791,
	     1.104167962462229,
	     1.062935513387498,
	     1.014676036983155,
	     0.959290764702396,
	     0.896687859619188,
	     0.826782162849412,
	     0.749494944197282,
	     0.664753657212279,
	     0.572491698817521,
	     0.472648173647942,
	     0.365167663215605};
    }else
      return{};
  }

};

    }}}
#endif
