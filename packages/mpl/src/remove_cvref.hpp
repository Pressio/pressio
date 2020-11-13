
#ifndef MPL_VARIADIC_REMOVE_CVREF_HPP_
#define MPL_VARIADIC_REMOVE_CVREF_HPP_

namespace pressio{ namespace mpl{

template<class T>
struct remove_cvref;

template<class T>
struct remove_cvref{
  using type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;
};

template<class T>
using remove_cvref_t = typename remove_cvref<T>::type;

}} // namespace pressio::mpl

#endif  // MPL_VARIADIC_REMOVE_CVREF_HPP_
