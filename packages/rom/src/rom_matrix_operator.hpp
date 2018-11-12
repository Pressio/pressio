
#ifndef ROM_MATRIX_OPERATOR_HPP_
#define ROM_MATRIX_OPERATOR_HPP_

#include "rom_operator_base.hpp"
#include "../../CORE_OPS"

namespace rompp{ namespace rom{

template<typename operator_type,
	 core::meta::enable_if_t<
	   core::meta::is_core_matrix_wrapper<
	     operator_type>::value
	   > * = nullptr
	 >
class MatrixOperator
  : public OperatorBase<MatrixOperator<operator_type>>{
      
public:
  MatrixOperator() = delete;
  explicit MatrixOperator(operator_type & opIn)
    : op_(&opIn){}
  ~MatrixOperator() = default;

private:
  friend OperatorBase<MatrixOperator<operator_type>>;
  operator_type * op_;

};//end class

}} // end namespace rompp::rom
#endif





  // //-------------------------------
  // //----      APPLY AS IS      ----
  // //-------------------------------
  // template <typename T, 
  //    core::meta::enable_if_t<
  //      core::details::traits<T>::wrapped_package_identifier ==
  //      core::details::WrappedPackageIdentifier::Eigen 
  //      > * = nullptr
  //    >
  // auto applyImpl(const T & X){
  //   return core::ops::product(*op_, X);
  // }
  // //---------------------------------

  // template <typename T1,
  // 	    typename T2, 
  //    core::meta::enable_if_t<
  //      (core::details::traits<T1>::wrapped_package_identifier ==
  // 	core::details::WrappedPackageIdentifier::Eigen ) and
  //      (core::details::traits<T2>::wrapped_package_identifier ==
  // 	core::details::WrappedPackageIdentifier::Eigen )
  //      > * = nullptr
  //    >
  // void applyImpl(const T1 & X, T2 & Y){
  //  core::ops::product(*op_, X, Y);
  // }
  // //---------------------------------
  
  // //---------------------------
  // //----     TRANSPOSE     ----
  // //---------------------------
  // template <typename T, 
  //    core::meta::enable_if_t<
  //      (core::details::traits<T>::wrapped_package_identifier ==
  // 	core::details::WrappedPackageIdentifier::Eigen ) 
  //      > * = nullptr
  //    >
  // auto applyTransposeImpl(const T & X){
  //   using opT_t = decltype(*op_);
  //   opT_t opT(*op_->data().transpose());
  //   return opT * (*X.data()); //core::ops::dot(*op_, X);
  // }
  // //---------------------------------

  // template <typename T1, typename T2,
  //    core::meta::enable_if_t<
  //      (core::details::traits<T1>::wrapped_package_identifier ==
  // 	core::details::WrappedPackageIdentifier::Eigen ) and
  //      (core::details::traits<T2>::wrapped_package_identifier ==
  // 	core::details::WrappedPackageIdentifier::Eigen )
  //      > * = nullptr
  //    >
  // void applyTransposeImpl(const T1 & X, T2 & Y){
  //   using opT_t = decltype(*op_);
  //   opT_t opT(*op_->data().transpose());
  //   *Y.data() = opT * (*X.data()); //core::ops::dot(*op_, X);
  //   //core::ops::dot(*op_, X, Y);
  // }
  // //---------------------------------
