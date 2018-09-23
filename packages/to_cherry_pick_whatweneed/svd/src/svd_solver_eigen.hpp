
#ifndef SVD_SOLVER_EIGEN_HPP_
#define SVD_SOLVER_EIGEN_HPP_

#include "svd_solver_generic_base.hpp"
#include "matrix/core_matrix_traits.hpp"

#include <Eigen/SVD>


namespace rompp{ 
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
  private:
    using derived_t = solver<matrix_type, which_impl>;
    using sc_t = typename svd::details::traits<derived_t>::scalar_t;
    using wrap_solver_t = typename svd::details::traits<derived_t>::wrapped_solver_t;
    using native_matrix_t = typename svd::details::traits<derived_t>::native_matrix_t;
    using u_matrix_type = typename svd::details::traits<derived_t>::u_matrix_t;
    using v_matrix_type = typename svd::details::traits<derived_t>::v_matrix_t;
   
  private:
    wrap_solver_t svdObj_;
    u_matrix_type U_;
    v_matrix_type V_;
    
  public:
    solver(){}
    ~solver(){}
    //-----------------------------------

    template <typename matrix_in_type>
    typename std::enable_if<std::is_same<matrix_in_type, matrix_type>::value>::type
    computeImpl(const matrix_in_type & mat){
      svdObj_.compute(*mat.view(), Eigen::ComputeThinU | Eigen::ComputeThinV );
      U_ = u_matrix_type( svdObj_.matrixU() );
      V_ = v_matrix_type( svdObj_.matrixV() );
    };

    const u_matrix_type & leftSingularVectorsImpl() const{
      return U_;
    };
    const v_matrix_type & rightSingularVectorsImpl() const{
      return V_;
    };
    
    // const  & singularValuesImpl(){
    // this needs to be fixed because we need to return values someway
    //   auto & singval = svdObj_.singularValues();
    //   std::cout << "sval:" << std::endl << singval << std::endl;
    //   //      return svdObj_.
    // };
    
    wrap_solver_t const & getConstRefToWrappedSolverImpl() const {
      return svdObj_;
    };
  
  };

  
}//end namespace svd  
}//end namespace rompp
#endif

