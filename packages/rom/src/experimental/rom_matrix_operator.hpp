
#ifndef ROM_MATRIX_OPERATOR_HPP_
#define ROM_MATRIX_OPERATOR_HPP_

#include "rom_operator_base.hpp"
#include "../../../CORE_OPS"

namespace rompp{ namespace rom{

template<typename operator_type,
	 core::meta::enable_if_t<
	   core::meta::is_core_matrix_wrapper<
	     operator_type>::value
	   > * = nullptr
	 >
class MatrixOperator
  : public OperatorBase<
      MatrixOperator<operator_type>>{

public:
  MatrixOperator() = delete;
  explicit MatrixOperator(operator_type & opIn)
    : op_(&opIn){}
  ~MatrixOperator() = default;

private:
  friend OperatorBase<MatrixOperator<operator_type> >;
  operator_type * op_;

};//end class

}} // end namespace rompp::rom
#endif
