
#ifndef ROM_MATRIX_OPERATOR_HPP_
#define ROM_MATRIX_OPERATOR_HPP_

#include "rom_fwd.hpp"
#include "../../ALGEBRA_OPS"

namespace rompp{ namespace rom{

template<typename wrapped_type>
class MatrixOperator<
  wrapped_type,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_algebra_matrix_wrapper<
      wrapped_type>::value
    >
  >{

public:
  MatrixOperator() = delete;

  explicit MatrixOperator(wrapped_type & opIn)
    : op_(&opIn){}

  ~MatrixOperator() = default;

private:
  wrapped_type * op_ = {};

};//end class

}} // end namespace rompp::rom
#endif
