
#ifndef ROM_MULTI_VECTOR_OPERATOR_HPP_
#define ROM_MULTI_VECTOR_OPERATOR_HPP_

#include "rom_fwd.hpp"
#include "../../CONTAINERS_OPS"

namespace pressio{ namespace rom{

template <typename wrapped_type>
class MultiVectorOperator <
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper<
      wrapped_type>::value
    >
  >{

  using this_t = MultiVectorOperator<wrapped_type>;
  const wrapped_type & op_ = {};

public:
  MultiVectorOperator() = delete;

  explicit MultiVectorOperator(const wrapped_type & opIn)
    : op_(opIn){}

  ~MultiVectorOperator() = default;

public:
  const wrapped_type * getPtrToOperator() const{
    return &op_;
  }

  const wrapped_type & getRefToOperator() const{
    return op_;
  }

  /* Y = X * op: return Y */
  template<typename T>
  auto applyRight(const T & X) const
    -> decltype(containers::ops::product(X, op_)){
    return containers::ops::product(X, op_);
  }

  /* Y = X * op: Y passed */
  template<typename T1, typename T2>
  void applyRight(const T1 & X, T2 & Y)  const{
    containers::ops::product(X, op_, Y);
  }


  /* Y = op_ * X : return Y */
  template <typename T>
  auto apply(const T & X) const
    -> decltype(containers::ops::product( std::declval<wrapped_type>(),
				    std::declval<T>() )){
    return containers::ops::product(op_, X);
  }

  /* Y = op_ * X : pass Y */
  template <typename T1,
	    typename T2,
     ::pressio::mpl::enable_if_t<
       containers::meta::is_vector_wrapper<T1>::value &&
       containers::meta::is_vector_wrapper<T2>::value
       > * = nullptr
     >
  void apply(const T1 & X, T2 & Y) const{
    // op_: multivector of size m,n
    // X: vector of size n,1
    // Y: vector of size m,1
    containers::ops::product(op_, X, Y);
  }


  /* Y = op^T * X: return Y */
  template <typename T,
     ::pressio::mpl::enable_if_t<
       containers::meta::is_vector_wrapper<T>::value or
       containers::meta::is_multi_vector_wrapper<T>::value
       > * = nullptr
     >
  auto applyTranspose(const T & X) const
    -> decltype(containers::ops::dot( std::declval<wrapped_type>(),
				std::declval<T>() )){
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    return containers::ops::dot(op_, X);
  }

  /* Y = op^T * X: pass Y */
  template <typename T1, typename T2,
     ::pressio::mpl::enable_if_t<
       containers::meta::is_vector_wrapper<T1>::value or
       containers::meta::is_multi_vector_wrapper<T1>::value
       > * = nullptr
     >
  void applyTranspose(const T1 & X, T2 & Y) const{
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    containers::ops::dot(op_, X, Y);
  }

};//end class

}} // end namespace pressio::rom
#endif
