
#ifndef ROM_IS_LEGITIMATE_MODEL_FOR_GALERKIN_HPP_
#define ROM_IS_LEGITIMATE_MODEL_FOR_GALERKIN_HPP_

#include "../../../ode/src/meta/ode_is_legitimate_model_for_explicit_ode.hpp"

namespace pressio{ namespace rom{ namespace meta {

template< typename T, typename enable = void >
struct is_legitimate_model_for_galerkin : std::false_type{};

template<typename T>
struct is_legitimate_model_for_galerkin<
  T,
  mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_model_for_explicit_ode<T>::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::meta
#endif
