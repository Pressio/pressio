
#ifndef ROM_MULTI_VECTOR_OPERATOR_HPP_
#define ROM_MULTI_VECTOR_OPERATOR_HPP_

//#include "rom_operator_base.hpp"
#include "rom_forward_declarations.hpp"
#include "../../CORE_ALL"

namespace rompp{ namespace rom{

template <typename operator_type>
class MultiVectorOperator < operator_type, 
			    core::meta::enable_if_t<
			      core::meta::is_core_multi_vector_wrapper<
				operator_type>::value
			      >
			    >
{

public:
  
  MultiVectorOperator() = delete;

  explicit MultiVectorOperator(const operator_type & opIn)
    : op_(&opIn){}

  ~MultiVectorOperator() = default;

private:
    using this_t = MultiVectorOperator<operator_type>;
  const operator_type * op_ = {};

public:

  const operator_type * getOperator() const{
    return op_;
  }

  //-------------------------------
  //----     APPLY RIGHT      -----
  //-------------------------------
  template<typename T>
  auto applyRight(const T & X)
    -> decltype(core::ops::product(X, *op_)){
    return core::ops::product(X, *op_);
  }

  template<typename T1, typename T2>
  void applyRight(const T1 & X, T2 & Y){
    core::ops::product(X, *op_, Y);
  }

  template <typename T>
  auto gigi(const T & X) -> double{
    auto gg = core::ops::product(*op_, X);
    return 0.0;
  }


  //-------------------------------
  //----     APPLY AS IS      -----
  //-------------------------------
  template <typename T
    /*core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value
       > * = nullptr */
     >
  auto apply(const T & X)
    -> decltype(core::ops::product( std::declval<operator_type>(),
				    std::declval<T>() )){
    return core::ops::product(*op_, X);
  }

  template <typename T1,
	    typename T2,
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T1>::value &&
       core::meta::is_core_vector_wrapper<T2>::value
       > * = nullptr
     >
  void apply(const T1 & X, T2 & Y){
    // op_: multivector of size m,n
    // X: vector of size n,1
    // Y: vector of size m,1
    core::ops::product(*op_, X, Y);
  }

  //---------------------------
  //----     TRANSPOSE     ----
  //---------------------------
  template <typename T,
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T>::value or
       core::meta::is_core_multi_vector_wrapper<T>::value
       > * = nullptr
     >
  auto applyTranspose(const T & X)
    -> decltype(core::ops::dot( std::declval<operator_type>(),
				std::declval<T>() )){
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    return core::ops::dot(*op_, X);
  }
  //---------------------------------

  template <typename T1, typename T2,
     core::meta::enable_if_t<
       core::meta::is_core_vector_wrapper<T1>::value or
       core::meta::is_core_multi_vector_wrapper<T1>::value
       > * = nullptr
     >
  void applyTranspose(const T1 & X, T2 & Y){
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    core::ops::dot(*op_, X, Y);
  }

};//end class

}} // end namespace rompp::rom
#endif
