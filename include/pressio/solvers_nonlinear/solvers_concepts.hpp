/*
//@HEADER
// ************************************************************************
//
// solvers_admissible_linear_solver_for_nonlinear_least_squares.hpp
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

#ifndef SOLVERS_NONLINEAR_CONCEPTS_HPP_
#define SOLVERS_NONLINEAR_CONCEPTS_HPP_

namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept SystemWithResidual =
     std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::residual_type>
  && Traits<typename T::state_type>::rank == 1
  && Traits<typename T::residual_type>::rank == 1
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::residual_type & r)
  {
    { A.createState()       } -> std::same_as<typename T::state_type>;
    { A.createResidual()    } -> std::same_as<typename T::residual_type>;
    { A.residual(state, r)  } -> std::same_as<void>;
  };

template <class T>
concept SystemWithResidualAndJacobian =
     SystemWithResidual<T>
  && std::copy_constructible<typename T::jacobian_type>
  && Traits<typename T::jacobian_type>::rank == 2
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::jacobian_type & j)
  {
    { A.createJacobian()    } -> std::same_as<typename T::jacobian_type>;
    { A.jacobian(state, j)  } -> std::same_as<void>;
  };

template <class T>
concept SystemWithFusedResidualAndJacobian =
     SystemWithResidual<T>
  && std::copy_constructible<typename T::jacobian_type>
  && Traits<typename T::jacobian_type>::rank == 2
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::residual_type & r,
	      typename T::jacobian_type & j,
	      bool computeJacobian)
  {
    { A.createJacobian()    } -> std::same_as<typename T::jacobian_type>;
    { A.residualAndJacobian(state, r, j)  } -> std::same_as<void>;
  };

template <class T>
concept RealValuedSystemWithResidual =
  SystemWithResidual<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::residual_type> >;

template <class T>
concept RealValuedSystemWithResidualAndJacobian =
     RealValuedSystemWithResidual<T>
  && SystemWithResidualAndJacobian<T>
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >;

template <class T>
concept RealValuedSystemWithFusedResidualAndJacobian =
     RealValuedSystemWithResidual<T>
  && SystemWithFusedResidualAndJacobian<T>
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >;


#else

template<class T, class enable = void>
struct SystemWithResidual : std::false_type{};

template<class T>
struct SystemWithResidual<
  T,
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::residual_type>::value
    && Traits<typename T::state_type>::rank == 1
    && Traits<typename T::residual_type>::rank == 1
    //
    && ::pressio::nonlinearsolvers::has_const_create_state_method_return_result<
      T, typename T::state_type>::value
    && ::pressio::nonlinearsolvers::has_const_create_residual_method_return_result<
      T, typename T::residual_type>::value
    //
    && ::pressio::nonlinearsolvers::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type>::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct SystemWithResidualAndJacobian : std::false_type{};

template<class T>
struct SystemWithResidualAndJacobian<
  T,
  mpl::enable_if_t<
    SystemWithResidual<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && Traits<typename T::jacobian_type>::rank == 2
    //
    && ::pressio::nonlinearsolvers::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type>::value
    && ::pressio::nonlinearsolvers::has_const_jacobian_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::jacobian_type>::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct SystemWithFusedResidualAndJacobian : std::false_type{};

template<class T>
struct SystemWithFusedResidualAndJacobian<
  T,
  mpl::enable_if_t<
    SystemWithResidual<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && Traits<typename T::jacobian_type>::rank == 2
    //
    && ::pressio::nonlinearsolvers::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type>::value
    && ::pressio::nonlinearsolvers::has_const_residualandjacobian_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type, typename T::jacobian_type>::value
   >
  > : std::true_type{};


template<class T, class = void> struct RealValuedSystemWithResidual : std::false_type{};
template<class T> struct RealValuedSystemWithResidual<
  T,
  mpl::enable_if_t<
    SystemWithResidual<T>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::residual_type> >::value
    >
  > : std::true_type{};

template<class T, class = void> struct RealValuedSystemWithResidualAndJacobian : std::false_type{};
template<class T> struct RealValuedSystemWithResidualAndJacobian<
  T,
  mpl::enable_if_t<
       RealValuedSystemWithResidual<T>::value
    && SystemWithResidualAndJacobian<T>::value
    && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
    >
  > : std::true_type{};


template<class T, class = void> struct RealValuedSystemWithFusedResidualAndJacobian : std::false_type{};
template<class T> struct RealValuedSystemWithFusedResidualAndJacobian<
  T,
  mpl::enable_if_t<
       RealValuedSystemWithResidual<T>::value
    && SystemWithFusedResidualAndJacobian<T>::value
    && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
    >
  > : std::true_type{};

#endif //PRESSIO_ENABLE_CXX20

//
// auxiliary stuff
//
#ifdef PRESSIO_ENABLE_CXX20
template <class T>
requires (RealValuedSystemWithResidual<T>
       || RealValuedSystemWithResidualAndJacobian<T>
       || RealValuedSystemWithFusedResidualAndJacobian<T>)
struct scalar_of{
  using type = scalar_trait_t< typename T::state_type >;
};
#else
template <class T, class = void> struct scalar_of;
template <class T>
struct scalar_of<
  T, mpl::enable_if_t<
       RealValuedSystemWithResidual<T>::value
       || RealValuedSystemWithResidualAndJacobian<T>::value
       || RealValuedSystemWithFusedResidualAndJacobian<T>::value
       > >
{
  using type = scalar_trait_t< typename T::state_type >;
};
#endif

template <class T> using scalar_of_t = typename scalar_of<T>::type;

}} // end namespace pressio::nonlinearsolvers

#endif  // SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_
