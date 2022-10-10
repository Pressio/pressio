/*
//@HEADER
// ************************************************************************
//
// ode_steppable.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ODE_ADVANCERS_CONSTRAINTS_ODE_STEPPABLE_WITH_ARGS_HPP_
#define ODE_ADVANCERS_CONSTRAINTS_ODE_STEPPABLE_WITH_ARGS_HPP_

#ifdef PRESSIO_ENABLE_CXX20

// this is here to use in the doc
//via rst literal include directive
namespace pressio{ namespace ode{

template <class T, class AuxT, class ...Args>
concept SteppableWithAuxiliaryArgs =
  /* must have nested aliases */
  requires(){
    typename T::independent_variable_type;
    typename T::state_type;
  }
  && requires(T & A,
	      typename T::state_type & state,
	      const ::pressio::ode::StepStartAt<typename T::independent_variable_type> & startAt,
	      const ::pressio::ode::StepCount & stepNumber,
	      const ::pressio::ode::StepSize<typename T::independent_variable_type> & dt,
	      AuxT && aux,
	      Args && ... args)
  {
    A(state, startAt, stepNumber, dt,
      std::forward<AuxT>(aux), std::forward<Args>(args)...);
  };

}} // end namespace pressio::ode

namespace pressio{ namespace ode{

template <class T, class AuxT, class ...Args>
concept StronglySteppableWithAuxiliaryArgs =
  SteppableWithAuxiliaryArgs<T, AuxT, Args...>;

}} // end namespace pressio::ode


/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else

namespace pressio{ namespace ode{ namespace impl{

/*
  for any stepper we need to ensure its operator() accepts the
  state by **non-const** reference. One (maybe only) way to detect
  that is check that if operator() can bind an rvalue state.
  Basically if we tried to bind an rvalue to an lvalue reference,
  it should fail, and if it fails that is good because it is what we want.

  In brief, we want the stepper to have this:
  class{
      void operator()(state_type & state, ...)
  };

  NOT this:
  class{
      void operator()(state_type state, ...)
  };
  or:
  class{
      void operator()(state_type && state, ...)
  };
 */

template <class T, class AuxT, class ...Args>
struct variadic_stepper_accepting_lvalue_state : std::true_type{};

template <class T, class AuxT, class ...Args>
struct variadic_stepper_accepting_lvalue_state<
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval<AuxT>(), std::declval<Args>()...
	)
       )
      >::value
    >, T, AuxT, Args...
  > : std::false_type{};

} //end impl

template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs : std::false_type{};

template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs<
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_independent_variable_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval<typename T::state_type & >(),
	std::declval<::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval<::pressio::ode::StepCount >(),
	std::declval<::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval< AuxT >(), std::declval<Args>()...
	)
       )
      >::value
    && impl::variadic_stepper_accepting_lvalue_state<void, T, AuxT, Args...>::value
    >,
  T, AuxT, Args...
  > : std::true_type{};


template <class T, class AuxT, class ...Args>
using StronglySteppableWithAuxiliaryArgs = SteppableWithAuxiliaryArgs<T, AuxT, Args...>;

}} // end namespace pressio::ode
#endif //end PRESSIO_ENABLE_CXX20

#endif  // ODE_ADVANCERS_CONSTRAINTS_ODE_STEPPABLE_WITH_ARGS_HPP_
