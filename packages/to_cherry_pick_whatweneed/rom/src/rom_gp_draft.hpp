
#ifndef ROM_GP_DRAFT_HPP_
#define ROM_GP_DRAFT_HPP_

#include "rom_ConfigDefs.hpp"
#include "svd_solver_traits.hpp"

namespace rom{

template <typename state_type,
	  typename app_type,
	  typename svdsolve_type>
class GP
{
private:
  using lsv_matrix_type = typename svd::details::traits<svdsolve_type>::u_matrix_t;    
  
  state_type * myY_;
  app_type * appPtr_;
  svdsolve_type * svd_;
  state_type Udoty_;
  lsv_matrix_type lsv_;
  lsv_matrix_type lsvT_;
  state_type dudtRed_;

public:
  GP(state_type & src,
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
    std::cout << "GP y0 " << std::endl;
    for (int i=0; i<Udoty_.size(); i++)
      std::cout << (*myY_)[i] << " ";
    std::cout << std::endl;
    std::cout << "--------- " << std::endl;

  // for (int i=0; i<Udoty_.size(); i++)
    //   std::cout << Udoty_[i] << std::endl;
    dudtRed_.resize( lsv_.rows() );
  }
  ~GP(){}

  
  void operator() ( const state_type & u,
		    state_type & dudt,
		    const double t )
  {
    // V * u
    u.matMultiply(lsv_, Udoty_);
    dudtRed_.resize( lsv_.rows() );
    (*appPtr_)( *Udoty_.view(), dudtRed_.getNonConstRefToData(), t);
    dudtRed_.matMultiply(lsvT_, dudt);
  };

};

}//end namespace rom

#endif

