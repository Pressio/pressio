
#ifndef SVD_SOLVER_EIGEN_HPP_
#define SVD_SOLVER_EIGEN_HPP_

#include "svd_solver_generic_base.hpp"
#include "matrix/core_matrix_traits.hpp"

#include <Eigen/SVD>


namespace svd{

     
  template <typename matrix_type,
	    svdKind which_impl>
  class solver<matrix_type, which_impl, 
	       typename std::enable_if< core::details::traits<matrix_type>::isEigen==1 &&
					(which_impl == svdKind::EigenJacobi ||
					 which_impl == svdKind::EigenBDCSVD)
					>::type
	       >
    : public solverGenericBase< solver<matrix_type, which_impl> >
  {
  public:
    using derived_t = solver<matrix_type, which_impl>;
    using sc_t = typename svd::details::traits<derived_t>::scalar_t;
    using wrap_solver_t = typename svd::details::traits<derived_t>::wrapped_solver_t;
        
  private:
    wrap_solver_t svdObj_;
      
  public:
    solver(){}
    ~solver(){}
    //-----------------------------------

    wrap_solver_t const & getConstRefToWrappedSolverImpl() const {
      return svdObj_;
    };
  
  };


  
}//end namespace svd
  
#endif

