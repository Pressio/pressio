
#ifndef rom_admissible_system_velocity_api_for_galerkin_HPP_
#define rom_admissible_system_velocity_api_for_galerkin_HPP_

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename enable = void>
struct rom_admissible_system_velocity_api_for_galerkin : std::false_type{};

template<typename T>
struct rom_admissible_system_velocity_api_for_galerkin<
  T,
  mpl::enable_if_t<  	
    ::pressio::ode::meta::admissible_system_explicit_ode<T>::value
    >
  > : std::true_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<>
// struct rom_admissible_system_velocity_api_for_galerkin<pybind11::object>
//   : std::true_type{};
// #endif

}}} // namespace pressio::rom::meta
#endif
