
#ifndef CORE_MATRIX_BASE_MATRIX_GENERIC_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_GENERIC_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixGenericBase
  : private core::details::CrtpBase<MatrixGenericBase<derived_type>>
{
private:
  using wrap_t = typename details::traits<derived_type>::wrapped_t;

public:
  wrap_t const * data() const {
    return this->underlying().dataImpl();
  };

  wrap_t * data() {
    return this->underlying().dataImpl();
  };
  
private:
  friend derived_type;
  friend core::details::CrtpBase<MatrixGenericBase<derived_type>>;
  MatrixGenericBase() = default;
  ~MatrixGenericBase() = default;
  
};//end class  
} // end namespace core
#endif
