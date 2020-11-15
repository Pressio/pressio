
#ifndef MPL_REMOVE_REFERENCE_HPP_
#define MPL_REMOVE_REFERENCE_HPP_

namespace pressio{ namespace mpl{

template<class T>
struct remove_reference;

template<class T>
struct remove_reference{
  using type = typename std::remove_reference<T>::type;
};

template<class T>
using remove_reference_t = typename remove_reference<T>::type;

}} // namespace pressio::mpl

#endif  // MPL_REMOVE_REFERENCE_HPP_
