
#ifndef CORE_MULTI_VECTOR_OPERATOR_HPP_
#define CORE_MULTI_VECTOR_OPERATOR_HPP_

#include "./core_operator_base.hpp"
#include "../multi_vector/core_multi_vector_traits.hpp"

namespace core{

template<typename operator_type,
	 core::meta::enable_if_t<
	   core::meta::is_core_multi_vector_wrapper<
	     operator_type>::value
	   > * = nullptr
	 >
class MultiVectorOperator
  : public OperatorBase<
      MultiVectorOperator<operator_type>>{

  template <typename operand_type,
	    typename result_type>
  void applyImpl(const operand_type & X,
		 result_type & Y,
		 bool useTranspose = false)
  {
    std::cout << "MultiVectorOperator" << std::endl;
    // do something
  }

public:
  MultiVectorOperator() = delete;

  explicit MultiVectorOperator(operator_type & opIn)
    : op_(&opIn){}

  ~MultiVectorOperator() = default;

private:
  friend OperatorBase<MultiVectorOperator<operator_type> >;

  operator_type * op_;

};//end class

} // end namespace core
#endif
