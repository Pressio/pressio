
#ifndef ROMPP_MPL_IDENTITY_HPP
#define ROMPP_MPL_IDENTITY_HPP

namespace rompp{ namespace mpl {

/**
 * \class identity
 * \brief Returns the argument passed
 */
template<class T> struct identity{
  using type = T;
};

}} // namespace rompp::mpl
#endif 
