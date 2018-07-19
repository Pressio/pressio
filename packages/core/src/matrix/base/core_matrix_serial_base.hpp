
#ifndef CORE_MATRIX_BASE_MATRIX_SERIAL_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SERIAL_BASE_HPP_

#include "../core_matrix_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class MatrixSerialBase
  : private core::details::CrtpBase<MatrixSerialBase<derived_type>>
{

private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;  

public:
  ord_t rows() const{
    return this->underlying().rowsImpl();
  }

  ord_t cols() const{
    return this->underlying().colsImpl();
  }

  void resize(ord_t nrows, ord_t ncols){
    this->underlying().resizeImpl(nrows, ncols);
  }

private:  
  friend derived_type;
  friend core::details::CrtpBase<MatrixSerialBase<derived_type>>;

  MatrixSerialBase() = default;
  ~MatrixSerialBase() = default;
  
};//end class
  
} // end namespace core
#endif
