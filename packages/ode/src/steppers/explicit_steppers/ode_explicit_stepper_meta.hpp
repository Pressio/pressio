
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPERS_META_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPERS_META_HPP_

#include "../../ode_basic_meta.hpp"
#include "../../policies/meta/ode_explicit_policies_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {
  
//---------------------------------------------------------------
// pair of types that matches:
//  state type | model type
//  model type  |  state type
template <typename T1, typename T2, typename enable = void>
struct is_expl_type_pair_state_and_model : std::false_type{};

template <typename T1, typename T2>
struct is_expl_type_pair_state_and_model<
  T1,T2,
  core::meta::enable_if_t< 
    ode::meta::is_legitimate_explicit_state_type<T1>::value and
    ode::meta::is_legitimate_model_for_explicit_ode<T2>::value
    > 
  > : std::true_type{
  using state_type = T1;
  using model_type = T2;
};

template <typename T1, typename T2>
struct is_expl_type_pair_state_and_model<
  T1,T2,
  core::meta::enable_if_t< 
    ode::meta::is_legitimate_model_for_explicit_ode<T1>::value and 
    ode::meta::is_legitimate_explicit_state_type<T2>::value
    >
  > : std::true_type{
  using state_type = T2;
  using model_type = T1;
};

  
//---------------------------------------------------------------
// pair of types that matches:
//  state type | residual policy type
//  residual policy type  |  state type

template <typename T1, typename T2, typename enable = void>
struct is_expl_type_pair_state_and_respol : std::false_type{};

template <typename T1, typename T2>
struct is_expl_type_pair_state_and_respol<
  T1,T2,
  core::meta::enable_if_t< 
    ode::meta::is_legitimate_explicit_state_type<T1>::value and
    ode::meta::is_legitimate_explicit_residual_policy<T2>::value
    > 
  > : std::true_type{
  using state_type = T1;
  using residual_policy_type = T2;
};

template <typename T1, typename T2>
struct is_expl_type_pair_state_and_respol<
  T1,T2,
  core::meta::enable_if_t< 
    ode::meta::is_legitimate_explicit_residual_policy<T1>::value and 
    ode::meta::is_legitimate_explicit_state_type<T2>::value
    >
  > : std::true_type{
  using state_type = T2;
  using residual_policy_type = T1;
};


//---------------------------------------------------------------
// sequence containing types: state, res pol, model 
template <typename T1, typename T2, typename T3, typename enable = void>
struct is_expl_type_set_state_respol_model : std::false_type{};

// T1,T2 = state/res pol  and T3 = model_t
template <typename T1, typename T2, typename T3>
struct is_expl_type_set_state_respol_model<
  T1,T2,T3,
  core::meta::enable_if_t<
    is_expl_type_pair_state_and_respol<T1,T2>::value and 
    ode::meta::is_legitimate_model_for_explicit_ode<T3>::value
    > 
  > : std::true_type{

  using pair_t = is_expl_type_pair_state_and_respol<T1,T2>;
  using state_type = typename pair_t::state_type;
  using residual_policy_type = typename pair_t::residual_policy_type;
  using model_type = T3;
};


// T1,T3 = state/res pol and T2 = model_t
template <typename T1, typename T2, typename T3>
struct is_expl_type_set_state_respol_model<
  T1,T2,T3,
  core::meta::enable_if_t<
    is_expl_type_pair_state_and_respol<T1,T3>::value and 
    ode::meta::is_legitimate_model_for_explicit_ode<T2>::value
    > 
  > : std::true_type{

  using pair_t = is_expl_type_pair_state_and_respol<T1,T3>;
  using state_type = typename pair_t::state_type;
  using residual_policy_type = typename pair_t::residual_policy_type;
  using model_type = T2;
};


// T2,T3 = state/res pol and T1 = model_t
template <typename T1, typename T2, typename T3>
struct is_expl_type_set_state_respol_model<
  T1,T2,T3,
  core::meta::enable_if_t<
    is_expl_type_pair_state_and_respol<T2,T3>::value and 
    ode::meta::is_legitimate_model_for_explicit_ode<T1>::value
    > 
  > : std::true_type{

  using pair_t = is_expl_type_pair_state_and_respol<T2,T3>;
  using state_type = typename pair_t::state_type;
  using residual_policy_type = typename pair_t::residual_policy_type;
  using model_type = T1;
};
      
      
}}} // end namespace rompp::ode::meta
#endif
