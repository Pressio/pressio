
#ifndef rom_compose_lspg_HPP_
#define rom_compose_lspg_HPP_

#include "./impl/unsteady/continuous_time_api/rom_lspg_unsteady_problem_continuous_time_api.hpp"
#include "./impl/unsteady/discrete_time_api/rom_lspg_unsteady_problem_discrete_time_api.hpp"

namespace pressio{ namespace rom{ namespace lspg{

namespace impl{

struct Default{};
struct Preconditioned{};
struct Masked{};


template<typename tag, typename ...Args>
struct compose{
  using type = void;
};

//********************
//***** UNSTEADY *****
//********************

// unsteady default lspg continuous time API
template<typename stepper_tag, typename system_type, typename lspg_state_t, typename ...Args>
struct compose<
::pressio::rom::lspg::impl::Default, 
mpl::enable_if_t< 
::pressio::rom::concepts::continuous_time_system<system_type>::value and 
(std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Euler>::value or
std::is_same< stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value)
>, 
stepper_tag, system_type, lspg_state_t, Args...>
{
  using type = ::pressio::rom::lspg::impl::unsteady::ProblemContinuousTimeApi<
            ::pressio::rom::lspg::impl::unsteady::DefaultProblemTraitsContinuousTimeApi, 
            stepper_tag, system_type, lspg_state_t, Args...>;
};

// unsteady default lspg discrete time api
template<typename stepper_tag, typename system_type, typename lspg_state_t, typename ...Args>
struct compose<
::pressio::rom::lspg::impl::Default, 
mpl::enable_if_t< 
::pressio::rom::concepts::discrete_time_system<system_type>::value and 
std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value 
>, 
stepper_tag, system_type, lspg_state_t, Args...>
{
  using type = ::pressio::rom::lspg::impl::unsteady::ProblemDiscreteTimeApi<
            ::pressio::rom::lspg::impl::unsteady::DefaultProblemTraitsDiscreteTimeApi, 
            stepper_tag, system_type, lspg_state_t, Args...>;
};

//********************
//***** STEADY *****
//********************
// unsteady default lspg continuous time API
template<typename system_type, typename lspg_state_t, typename ...Args>
struct compose<
::pressio::rom::lspg::impl::Default, 
mpl::enable_if_t< 
::pressio::rom::concepts::steady_system<system_type>::value 
>, 
system_type, lspg_state_t, Args...>
{
  using type = ::pressio::rom::lspg::impl::unsteady::ProblemContinuousTimeApi<
            ::pressio::rom::lspg::impl::unsteady::DefaultProblemTraitsContinuousTimeApi, 
            stepper_tag, system_type, lspg_state_t, Args...>;
};

}// end impl namespace


template<typename ...Args>
using composeDefaultProblem = impl::compose<impl::Default, void, Args...>;

// template<typename ...Args>
// using composePreconditionedProblem = impl::compose<impl::Preconditioned, void, Args...>;

// template<typename ...Args>
// using composeMaskedProblem = impl::compose<impl::Masked, void, Args...>;

}}}
#endif
