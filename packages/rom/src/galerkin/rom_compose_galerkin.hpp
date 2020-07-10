
#ifndef rom_compose_galerkin_HPP_
#define rom_compose_galerkin_HPP_

#include "./impl/continuous_time_api/rom_galerkin_problem_continuous_time_api.hpp"
#include "./impl/discrete_time_api/rom_galerkin_problem_discrete_time_api.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

namespace impl{

struct Default{};

template<typename tag, typename ...Args>
struct compose{
  using type = void;
};

// default continuous time API
template<typename stepper_tag, typename system_type, typename galerkin_state_t, typename ...Args>
struct compose<
::pressio::rom::galerkin::impl::Default, 
mpl::enable_if_t< 
::pressio::rom::concepts::continuous_time_system<system_type>::value and 
(std::is_same< stepper_tag, ::pressio::ode::explicitmethods::Euler>::value or
std::is_same< stepper_tag, ::pressio::ode::explicitmethods::RungeKutta4>::value)
>, 
stepper_tag, system_type, galerkin_state_t, Args...>
{
  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
            ::pressio::rom::galerkin::impl::DefaultProblemTraitsContinuousTimeApi, 
            stepper_tag, system_type, galerkin_state_t, Args...>;
};

// default discrete time api
template<typename stepper_tag, typename system_type, typename galerkin_state_t, typename ...Args>
struct compose<
::pressio::rom::galerkin::impl::Default, 
mpl::enable_if_t< 
::pressio::rom::concepts::discrete_time_system<system_type>::value and 
std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value 
>, 
stepper_tag, system_type, galerkin_state_t, Args...>
{
  using type = ::pressio::rom::galerkin::impl::ProblemDiscreteTimeApi<
            ::pressio::rom::galerkin::impl::DefaultProblemTraitsDiscreteTimeApi, 
            stepper_tag, system_type, galerkin_state_t, Args...>;
};

}// end impl namespace

template<typename ...Args>
using composeDefaultProblem = impl::compose<impl::Default, void, Args...>;

}}}
#endif
