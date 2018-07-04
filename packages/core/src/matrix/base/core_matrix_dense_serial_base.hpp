
#ifndef CORE_MATRIX_DENSE_SERIAL_BASE_HPP_
#define CORE_MATRIX_DENSE_SERIAL_BASE_HPP_

#include "../core_matrix_traits.hpp"
#include "../../core_operators_base.hpp"

namespace core{
    
template<typename derived_type>
class matrixDenseSerialBase
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
  // the () operators below are placed here becasue they serve the
  // purpose of scripting operator for matrices, but this is just meant
  // to be for DENSE SERIAL matrices. So we leave them here.
  // Every dense serial matrix should define these.
  sc_t & operator() (ord_t row, ord_t col);
  sc_t const & operator() (ord_t row, ord_t col) const;

private:  
  friend derived_type;
   matrixDenseSerialBase() = default;
  ~matrixDenseSerialBase() = default;

private:
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
  
};//end class    
} // end namespace core
#endif
