

#ifndef rom_admissible_system_velocity_api_for_wls_HPP_
#define rom_admissible_system_velocity_api_for_wls_HPP_

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename enable = void>
struct rom_admissible_system_velocity_api_for_wls : std::false_type{};

template<typename T>
struct rom_admissible_system_velocity_api_for_wls<
  T,
  mpl::enable_if_t<
    ::pressio::rom::meta::admissible_system_velocity_api_unsteady_lspg<T>::value
    >
  > : std::true_type{};


}}} // namespace pressio::rom::meta
#endif
