
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
  
  // template <typename operand_type,
  // 	    typename result_type>
  // void applyImpl(const operand_type & X,
  // 		 result_type & Y,
  // 		 bool useTranspose = false){
  //   // do something
  // }
  // //---------------------------------
  
  template <typename T, 
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value 
       > * = nullptr
     >
  auto applyImpl(const T & X){

    //std::cout << "MultiVectorOperator" << std::endl;
    op_->data()->Print(std::cout);
    return core::ops::product(*op_, X);
  }
  //---------------------------------

  template <typename T, 
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value
       > * = nullptr
     >
  auto applyTranspImpl(const T & X){

    std::cout << "MultiVectorOperator trans" << std::endl;
    return core::ops::dot(*op_, X);
  }
  //---------------------------------
  
  
public:
  MultiVectorOperator() = delete;

  explicit MultiVectorOperator(operator_type & opIn)
    : op_(&opIn){}

  ~MultiVectorOperator() = default;

private:
  friend OperatorBase<MultiVectorOperator<operator_type> >;

  operator_type * op_;

};//end class

} // end namespace rom
}//end namespace rompp
#endif
