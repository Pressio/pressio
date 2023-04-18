
#ifndef ODE_CONCEPTS_SYSTEM_ALL_HPP_
#define ODE_CONCEPTS_SYSTEM_ALL_HPP_

#include "ode_has_const_discrete_residual_jacobian_method.hpp"

namespace pressio{ namespace ode{

template <class T>
concept OdeSystem =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::rhs_type & f)
  {
    { A.createState() } -> std::same_as<typename T::state_type>;
    { A.createRhs()   } -> std::same_as<typename T::rhs_type>;
    { A.rhs(state, evalValue, f) } -> std::same_as<void>;
  };

template <class T>
concept OdeSystemFusingRhsAndJacobian =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::jacobian_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::rhs_type & f,
	      std::optional<typename T::jacobian_type*> J)
  {
    { A.createState()    } -> std::same_as<typename T::state_type>;
    { A.createRhs()      } -> std::same_as<typename T::rhs_type>;
    { A.createJacobian() } -> std::same_as<typename T::jacobian_type>;
    { A.rhsAndJacobian(state, evalValue, f, J) } -> std::same_as<void>;
  };

template <class T>
concept OdeSystemFusingMassMatrixAndRhs =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::mass_matrix_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::mass_matrix_type & M,
	      typename T::rhs_type & f)
  {
    { A.createState()      } -> std::same_as<typename T::state_type>;
    { A.createRhs()        } -> std::same_as<typename T::rhs_type>;
    { A.createMassMatrix() } -> std::same_as<typename T::mass_matrix_type>;
    { A.massMatrixAndRhs(state, evalValue, M, f) } -> std::same_as<void>;
  };

template <class T>
concept CompleteOdeSystem =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::mass_matrix_type>
  && std::copy_constructible<typename T::jacobian_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::mass_matrix_type & M,
	      typename T::rhs_type & f,
	      std::optional<typename T::jacobian_type*> J)
  {
    { A.createState()      } -> std::same_as<typename T::state_type>;
    { A.createRhs()        } -> std::same_as<typename T::rhs_type>;
    { A.createMassMatrix() } -> std::same_as<typename T::mass_matrix_type>;
    { A.massMatrixAndRhsAndJacobian(state, evalValue, M, f, J) } -> std::same_as<void>;
  };

template <class T, int NumStates>
concept FullyDiscreteSystemWithJacobian =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::discrete_residual_type>
  && std::copy_constructible<typename T::discrete_jacobian_type>
  && requires(const T & A)
  {
    { A.createState()            } -> std::same_as<typename T::state_type>;
    { A.createDiscreteResidual() } -> std::same_as<typename T::discrete_residual_type>;
    { A.createDiscreteJacobian() } -> std::same_as<typename T::discrete_jacobian_type>;
  }
  /*todo: fix syntax */
  && ::pressio::ode::has_const_discrete_residual_jacobian_method<
    T, NumStates,
    typename ::pressio::ode::StepCount::value_type,
    typename T::independent_variable_type,
    typename T::state_type,
    typename T::discrete_residual_type,
    typename T::discrete_jacobian_type>::value;

//
// refine for real-valued case
//
template <class T>
concept RealValuedOdeSystem =
     OdeSystem<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T>
concept RealValuedOdeSystemFusingMassMatrixAndRhs =
  OdeSystemFusingMassMatrixAndRhs<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::mass_matrix_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T>
concept RealValuedOdeSystemFusingRhsAndJacobian =
     OdeSystemFusingRhsAndJacobian<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T>
concept RealValuedCompleteOdeSystem =
     CompleteOdeSystem<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::mass_matrix_type> >
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T, int NumStates>
concept RealValuedFullyDiscreteSystemWithJacobian =
     FullyDiscreteSystemWithJacobian<T, NumStates>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::discrete_residual_type> >
  && std::floating_point< scalar_trait_t<typename T::discrete_jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

//
// policy
//
template <class T>
concept ImplicitResidualJacobianPolicy =
  requires(){
    typename T::independent_variable_type;
    typename T::state_type;
    typename T::residual_type;
    typename T::jacobian_type;
  }
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
  && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
  && all_have_traits_and_same_scalar<
    typename T::state_type,
    typename T::residual_type,
    typename T::jacobian_type>::value
  && std::convertible_to<
    typename T::independent_variable_type,
    scalar_trait_t<typename T::state_type>>
  && requires(const T & A)
  {
    { A.createState()    } -> std::same_as<typename T::state_type>;
    { A.createResidual() } -> std::same_as<typename T::residual_type>;
    { A.createJacobian() } -> std::same_as<typename T::jacobian_type>;
  }
  && requires(const T & A,
	      StepScheme schemeName,
	      const typename T::state_type & predictedState,
	      const ImplicitStencilStatesDynamicContainer<typename T::state_type> & stencilStatesManager,
	      ImplicitStencilRightHandSideDynamicContainer<typename T::residual_type> & stencilVelocities,
	      const ::pressio::ode::StepEndAt<typename T::independent_variable_type> & rhsEvaluationTime,
	      ::pressio::ode::StepCount stepNumber,
	      const ::pressio::ode::StepSize<typename T::independent_variable_type> & dt,
	      typename T::residual_type & R,
	      std::optional<typename T::jacobian_type *> & J)
  {
    A(schemeName, predictedState, stencilStatesManager,
      stencilVelocities, rhsEvaluationTime, stepNumber, dt, R, J);
  };



//
// auxiliary stuff
//
template <class T, int n>
requires (   RealValuedOdeSystem<T>
	  || RealValuedOdeSystemFusingRhsAndJacobian<T>
	  || RealValuedOdeSystemFusingMassMatrixAndRhs<T>
	  || RealValuedCompleteOdeSystem<T>
	  || RealValuedFullyDiscreteSystemWithJacobian<T, n>
	 )
struct scalar_of{
  using type = scalar_trait_t< typename T::state_type >;
};

template <class T, int n=0> using scalar_of_t = typename scalar_of<T,n>::type;

}}
#endif
