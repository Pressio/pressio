
#ifndef ROMPP_MPL_ENABLE_IF_T_HPP_
#define ROMPP_MPL_ENABLE_IF_T_HPP_

#include <type_traits>

namespace rompp{ namespace mpl{

/*This allows allows for a shorter syntax:
    enable_if_t<a_condition, MyType>
  as opposed to:
    typename enable_if<a_certain_condition, MyType>::type
*/
template<bool condition, typename T = void>
using enable_if_t = typename std::enable_if<condition,T>::type;

}} // namespace rompp::mpl
#endif
