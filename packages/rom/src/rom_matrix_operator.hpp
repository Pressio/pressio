
#ifndef ROM_MATRIX_OPERATOR_HPP_
#define ROM_MATRIX_OPERATOR_HPP_

#include "rom_fwd.hpp"
#include "../../CONTAINERS_OPS"

namespace pressio{ namespace rom{

template<typename wrapped_type>
class MatrixOperator<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper<
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

}} // end namespace pressio::rom
#endif
