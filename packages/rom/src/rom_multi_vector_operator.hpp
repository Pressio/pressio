
#ifndef ROM_MULTI_VECTOR_OPERATOR_HPP_
#define ROM_MULTI_VECTOR_OPERATOR_HPP_

#include "rom_fwd.hpp"
#include "../../ALGEBRA_OPS"

namespace rompp{ namespace rom{

template <typename wrapped_type>
class MultiVectorOperator <
  wrapped_type,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_algebra_multi_vector_wrapper<
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
    -> decltype(algebra::ops::product(X, op_)){
    return algebra::ops::product(X, op_);
  }

  /* Y = X * op: Y passed */
  template<typename T1, typename T2>
  void applyRight(const T1 & X, T2 & Y)  const{
    algebra::ops::product(X, op_, Y);
  }


  /* Y = op_ * X : return Y */
  template <typename T>
  auto apply(const T & X) const
    -> decltype(algebra::ops::product( std::declval<wrapped_type>(),
				    std::declval<T>() )){
    return algebra::ops::product(op_, X);
  }

  /* Y = op_ * X : pass Y */
  template <typename T1,
	    typename T2,
     ::rompp::mpl::enable_if_t<
       algebra::meta::is_algebra_vector_wrapper<T1>::value &&
       algebra::meta::is_algebra_vector_wrapper<T2>::value
       > * = nullptr
     >
  void apply(const T1 & X, T2 & Y) const{
    // op_: multivector of size m,n
    // X: vector of size n,1
    // Y: vector of size m,1
    algebra::ops::product(op_, X, Y);
  }


  /* Y = op^T * X: return Y */
  template <typename T,
     ::rompp::mpl::enable_if_t<
       algebra::meta::is_algebra_vector_wrapper<T>::value or
       algebra::meta::is_algebra_multi_vector_wrapper<T>::value
       > * = nullptr
     >
  auto applyTranspose(const T & X) const
    -> decltype(algebra::ops::dot( std::declval<wrapped_type>(),
				std::declval<T>() )){
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    return algebra::ops::dot(op_, X);
  }

  /* Y = op^T * X: pass Y */
  template <typename T1, typename T2,
     ::rompp::mpl::enable_if_t<
       algebra::meta::is_algebra_vector_wrapper<T1>::value or
       algebra::meta::is_algebra_multi_vector_wrapper<T1>::value
       > * = nullptr
     >
  void applyTranspose(const T1 & X, T2 & Y) const{
    // multivector^T acts on vector = take dot of each row
    // op_^T: multivector of size n,m
    // X: vector of size m,1
    // Y: vector with results of all dots of size n,1
    algebra::ops::dot(op_, X, Y);
  }

};//end class

}} // end namespace rompp::rom
#endif
