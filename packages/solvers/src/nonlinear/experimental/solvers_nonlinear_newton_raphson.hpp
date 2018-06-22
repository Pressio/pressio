
#ifndef SOLVERS_NONLINEAR_NEWTON_RAPHSON_HPP_
#define SOLVERS_NONLINEAR_NEWTON_RAPHSON_HPP_

#include "solvers_ConfigDefs.hpp"
#include "linear/solvers_linear_eigen.hpp"

namespace solvers{
namespace experimental{
  

template <typename state_type,
	  typename residual_type,
	  typename jacobian_type,
	  typename linear_solver_type,
	  typename 
	  std::enable_if<core::details::traits<jacobian_type>::isMatrix==1 &&
			 core::details::traits<jacobian_type>::isEigen==1 &&
			 core::details::traits<residual_type>::isVector==1 &&
			 core::details::traits<residual_type>::isEigen==1 && 
			 core::details::traits<state_type>::isVector==1 &&
			 core::details::traits<state_type>::isEigen==1
			 >::type * = nullptr
	  >
class newtonRaphson
{
  using scalar_t = typename core::details::traits<state_type>::scalar_t;
  static constexpr int maxTrials = 50;
  static constexpr scalar_t tolf = 1e-12;
  static constexpr scalar_t neg1 = static_cast<scalar_t>(-1);
  static constexpr scalar_t zero = static_cast<scalar_t>(0);

public:
  newtonRaphson(const linear_solver_type & ls) : ls_(&ls){}
  ~newtonRaphson() = default;

  template<typename evaluator_type>
  void solve(state_type & y, evaluator_type & appObj) const
  {
    // create aux variables 
    state_type dy(y.size());
    state_type RE(y.size());
    jacobian_type JA(RE.size(), y.size());

    scalar_t errf = zero;
    scalar_t erry = zero;
    for (int step=0; step<maxTrials; step++)
    {
      //      std::cout << "STEP " << step << std::endl;
      // compute residual and jacobian for initial guess
      appObj.residual(y, RE);
      appObj.jacobian(y, JA);
      // std::cout << *RE.data() << std::endl;
      // std::cout << *JA.data() << std::endl;

      //check func convergence
      errf=zero;
      for (decltype(y.size()) i=0; i<y.size(); i++)
    	errf += std::abs(y[i]);
      //      std::cout << " erff " << errf << std::endl;
      if (errf < tolf)
    	break;

      // solve J dy = -F
      ls_->solve(JA, RE, dy);

      erry = zero;
      for (decltype(y.size()) i=0; i < y.size(); i++){
    	erry += std::abs(dy[i]);
    	y[i] -= dy[i];
      }
      //      std::cout << " erry " << erry << std::endl;
      if (erry < tolf)
    	break;
    };//end for
  };//end solve

private:
  linear_solver_type const * ls_;
  
};//end class
  
}//end namespace experimental
}//end namespace solvers
  
#endif
