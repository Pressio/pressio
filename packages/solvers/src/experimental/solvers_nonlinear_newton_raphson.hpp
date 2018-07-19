
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_NEWTON_RAPHSON_HPP_
#define SOLVERS_EXPERIMENTAL_NONLINEAR_NEWTON_RAPHSON_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_linear_eigen.hpp"

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
  static constexpr int maxTrials = 100;
  static constexpr scalar_t tolf = 1e-10;
  static constexpr scalar_t neg1 = static_cast<scalar_t>(-1);
  static constexpr scalar_t zero = static_cast<scalar_t>(0);

public:
  explicit newtonRaphson(linear_solver_type & ls)
    : ls_(&ls)
  {}
  ~newtonRaphson() = default;

  template<typename evaluator_type>
  void solve(state_type & y, evaluator_type & appObj)
  {
    // std::cout << "-------------" << std::endl;
    // std::cout << "NEWTONRAPHSON" << std::endl;
    // std::cout << "-------------" << std::endl;

    state_type dy_(y.size());
    state_type RE_(y.size());
    jacobian_type JA_(RE_.size(), y.size());
    
    scalar_t errf = zero;
    scalar_t erry = zero;
    for (int step=0; step<maxTrials; step++)
    {
      //      std::cout << "\nSTEP " << step << std::endl;
      // compute residual and jacobian for current guess      
      appObj.residual(y, RE_);
      appObj.jacobian(y, JA_);
      // std::cout << "RE" << std::endl;
      // for (int i=0; i<RE_.size();i++)
      //  	std::cout << std::setprecision(15) << RE_[i]  << " ";
      // std::cout<< "\n";
      //  std::cout << "JA" << std::endl;
      // std::cout << *JA_.data() << std::endl;
      // std::cout << "--------------------------" << std::endl;

      // //check func convergence
      // errf=zero;
      // for (decltype(y.size()) i=0; i<y.size(); i++)
      // 	errf += std::abs(y[i]);
      // //      std::cout << " erff " << errf << std::endl;
      // if (errf < tolf)
      // 	break;

      // solve J dy = F
      for (int i=0; i<dy_.size(); i++)
        dy_[i] = 0.0;
      ls_->solve(JA_, RE_, dy_);
      // std::cout << "dy" << std::endl;
      // for (int i=0; i<dy_.size();i++)
      //  	std::cout << std::setprecision(15) << dy_[i]  << " ";
      // std::cout<< "\n";

      erry = zero;
      for (decltype(y.size()) i=0; i < y.size(); i++){
    	erry += std::abs(dy_[i]);
    	y[i] -= dy_[i];
      }
      // std::cout << "new y" << std::endl;
      // for (int i=0; i<y.size();i++)
      //  	std::cout << std::setprecision(15) << y[i]  << " ";
      // std::cout<< "\n";
      if (erry < tolf)
      	break;
    };//end for

  }//end solve

private:
  linear_solver_type * ls_;
  // state_type dy_;
  // state_type RE_;
  // jacobian_type JA_;
  
};//end class
  
}//end namespace experimental
}//end namespace solvers
  
#endif
