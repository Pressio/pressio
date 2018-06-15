
#ifndef ROM_LSPG_MOCK_HPP_
#define ROM_LSPG_MOCK_HPP_

#include "rom_ConfigDefs.hpp"
#include "svd_solver_traits.hpp"

namespace rom{

template <typename state_type,
	  typename jacobian_type,
	  typename app_type,
	  typename svdsolve_type>
class lspg
{
private:
  using lsv_matrix_type = typename svd::details::traits<svdsolve_type>::u_matrix_t;
  state_type * myY_;
  app_type * appPtr_;
  svdsolve_type * svd_;
  state_type Udoty_;
  lsv_matrix_type lsv_;
  lsv_matrix_type lsvT_;
  // state_type dudtRed_;

public:
  lspg(state_type & src,
       app_type & obj,
       svdsolve_type & svd)
    : myY_(&src), appPtr_(&obj), svd_(&svd)
  {
    lsv_ = svd_->leftSingularVectors();
    lsv_.transpose(lsvT_);
    Udoty_.resize( lsvT_.cols() );
    // U^T * y0
    myY_->matMultiply(lsvT_, Udoty_);
    *myY_ = Udoty_;
    std::cout << "LSPG y0 " << std::endl;
    for (int i=0; i<Udoty_.size(); i++)
      std::cout << (*myY_)[i] << " ";
    std::cout << std::endl;
    std::cout << "--------- " << std::endl;

  }

  ~lspg(){}
  
  void operator() ( const state_type & u,
		    state_type & dudt,
		    jacobian_type & jac,
		    const double t )
  {
    
    // V * u
    u.matMultiply(lsv_, Udoty_);
    dudt.resize( lsv_.rows() );
    jac.getNonConstRefToData() = core::details::traits<jacobian_type>::wrapped_t::Zero(lsv_.rows(), lsv_.rows());
    
    (*appPtr_)( *Udoty_.view(), dudt.getNonConstRefToData(), jac.getNonConstRefToData(), t);    
  };
  
  void rescaleState(const state_type & uin,
		    state_type & uout){
    uout.resize( lsv_.rows() );
    uin.matMultiply(lsv_, uout);
  }
  
  void rescaleJacobian(const jacobian_type & jacfullin,
		       jacobian_type & JdotV){
    JdotV.getNonConstRefToData() = core::details::traits<jacobian_type>::wrapped_t::Zero(lsv_.rows(), lsv_.cols());
    //JdotV.resize( lsv_.rows(), lsv_.cols() );
    JdotV.getNonConstRefToData() = *jacfullin.view() * *lsv_.view();
  }
    
  
};

}//end namespace rom
#endif

