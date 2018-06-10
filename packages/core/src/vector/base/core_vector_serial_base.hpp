
#ifndef CORE_VECTOR_SERIAL_BASE_HPP_
#define CORE_VECTOR_SERIAL_BASE_HPP_

#include "../core_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class vectorSerialBase
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;
  
public:
  size_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(size_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
    
private:
  //friend class derived_type;
  vectorSerialBase(){}
  ~vectorSerialBase(){}
 
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };
  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };
   
};//end class
    
} // end namespace core
#endif


// template <typename wrapped_matrix_type, typename U>
// typename std::enable_if<std::is_same<U,der_t>::value>::type
// matMultiply(const matrix<wrapped_matrix_type> & matin,
// 		   U & result) const{
//   this->underlying().matMultiplyImpl(matin, result);
// };
