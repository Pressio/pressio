
#ifndef PRESSIO_MPL_IDENTITY_HPP
#define PRESSIO_MPL_IDENTITY_HPP

namespace pressio{ namespace mpl {

/**
 * \brief Returns the argument passed
 */
template<class T> struct identity{
  using type = T;
};

}} // namespace pressio::mpl
#endif 
