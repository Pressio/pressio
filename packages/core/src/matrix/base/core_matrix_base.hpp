
#ifndef CORE_MATRIX_BASE_MATRIX_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{
namespace core{
    
template<typename derived_type>
class MatrixBase
  : private core::details::CrtpBase<
  MatrixBase<derived_type>>{

  using traits_t = details::traits<derived_type>;
  using sc_t = typename traits_t::scalar_t;

public:
  void setIdentity(){
    this->underlying().setIdentityImpl();}

  bool isLowerTriangular(){
    return this->underlying().isLowerTriangularImpl();}
  
  bool isUpperTriangular(){
    return this->underlying().isUpperTriangularImpl();}
    
private:  
  friend derived_type;
  friend core::details::CrtpBase<MatrixBase<derived_type>>;
  MatrixBase() = default;
  ~MatrixBase() = default; 
};//end class
  
} // end namespace core
}//end namespace rompp
#endif
