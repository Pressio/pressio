
#ifndef CORE_MATRIX_BASE_MATRIX_GENERIC_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_GENERIC_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixGenericBase
  : private core::details::CrtpBase<MatrixGenericBase<derived_type>>
{

private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

public:
  wrap_t const * data() const {
    return this->underlying().dataImpl();
  };

  wrap_t * data() {
    return this->underlying().dataImpl();
  };

  void addToDiagonal(sc_t value) {
    return this->underlying().addToDiagonalImpl(value);
  };

  
private:
  friend derived_type;
  friend core::details::CrtpBase<MatrixGenericBase<derived_type>>;

  MatrixGenericBase() = default;
  ~MatrixGenericBase() = default;
  
};//end class
  
} // end namespace core
#endif
