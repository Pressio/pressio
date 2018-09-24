
#ifndef ROM_MULTI_VECTOR_OPERATOR_HPP_
#define ROM_MULTI_VECTOR_OPERATOR_HPP_

#include "rom_operator_base.hpp"
#include "../../../CORE_ALL"

namespace rompp{
namespace rom{

template<typename operator_type,
	 core::meta::enable_if_t<
	   core::meta::is_core_multi_vector_wrapper<
	     operator_type>::value
	   > * = nullptr
	 >
class MultiVectorOperator
  : public OperatorBase<
      MultiVectorOperator<operator_type>>{
  
  //---------------------------------

  template <typename T, 
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value 
       > * = nullptr
     >
  auto applyImpl(const T & X){
    return core::ops::product(*op_, X);
  }
  //---------------------------------

  template <typename T1,
	    typename T2, 
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T1>::value &&
       core::meta::is_core_vector_wrapper<T2>::value 
       > * = nullptr
     >
  void applyImpl(const T1 & X, const T2 & Y){
    core::ops::product(*op_, X);
  }

  //---------------------------------
  //---------------------------------
  //----      TRANSPOSE
  //---------------------------------
  //---------------------------------
  
  template <typename T, 
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value
       > * = nullptr
     >
  auto applyTranspImpl(const T & X){
    return core::ops::dot(*op_, X);
  }
  //---------------------------------
  
  
public:
  MultiVectorOperator() = delete;

  explicit MultiVectorOperator(const operator_type & opIn)
    : op_(&opIn){}

  ~MultiVectorOperator() = default;

private:
  friend OperatorBase<MultiVectorOperator<operator_type> >;

  const operator_type * op_;

};//end class

} // end namespace rom
}//end namespace rompp
#endif
