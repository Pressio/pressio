
#ifndef rom_continuous_time_system_hpp_
#define rom_continuous_time_system_hpp_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system : std::false_type{};

template<typename T>
struct continuous_time_system<
  T,
  mpl::enable_if_t<
    ::pressio::rom::concepts::continuous_time_explicit_system<T>::value or
    ::pressio::rom::concepts::continuous_time_implicit_system<T>::value 
    >
  > : std::true_type{};

}}} // namespace pressio::rom::concepts
#endif
