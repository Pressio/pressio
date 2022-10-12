/*
//@HEADER
// ************************************************************************
//
// rom_fom_system_continuous_time.hpp
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

#ifndef ROM_CONCEPTS_POSSIBLY_AFFINE_TRIAL_COLUMN_SUBSPACE_HPP_
#define ROM_CONCEPTS_POSSIBLY_AFFINE_TRIAL_COLUMN_SUBSPACE_HPP_

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept PossiblyAffineTrialColumnSubspace =
  /* must have nested aliases */
  requires(){
    typename T::reduced_state_type;
    typename T::basis_matrix_type;
    typename T::full_state_type;
  }
  /*
    requirements on the nested aliases
  */
  && pressio::is_vector_eigen<typename T::reduced_state_type>::value
  && std::copy_constructible< typename T::full_state_type>
  && std::copy_constructible< typename T::basis_matrix_type>
  && all_have_traits_and_same_scalar<
	typename T::reduced_state_type,
	typename T::full_state_type,
	typename T::basis_matrix_type>::value
  && Traits<typename T::full_state_type>::rank == 1
  && Traits<typename T::basis_matrix_type>::rank == 2
  /**/
  && std::copy_constructible<T>
  && !std::assignable_from<T&, T>
  && !std::assignable_from<T&, T&>
  /*
    compound requirements on the "create*" methods
  */
  && requires(const T & A,
	      const typename T::reduced_state_type & reducedState,
	      typename T::full_state_type & fullState)
  {
    {A.createReducedState() } -> std::same_as<typename T::reduced_state_type>;
    {A.createFullState()    } -> std::same_as<typename T::full_state_type>;

    {A.createFullStateFromReducedState(reducedState) }
	   -> std::same_as<typename T::full_state_type>;
  }
  /*
    other compound requirements
  */
  && requires(const T & A,
	      const typename T::reduced_state_type & reducedState,
	      typename T::full_state_type & fullState)
  {
    {A.mapFromReducedState(reducedState,
			   fullState) } -> std::same_as<void>;

    {A.basisOfTranslatedSpace()} -> std::same_as<const typename T::basis_matrix_type &>;
    {A.translationVector()     } -> std::same_as<const typename T::full_state_type &>;

    {A.basis()        } -> std::same_as<const typename T::basis_matrix_type &>;
    {A.dimension()    } -> std::integral;
    {A.isColumnSpace()} -> std::convertible_to<bool>;
    {A.isRowSpace()   } -> std::convertible_to<bool>;
  };
#endif // PRESSIO_ENABLE_CXX20

}} //end namespace pressio::rom



#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace rom{

template<class T, class enable = void>
struct PossiblyAffineTrialColumnSubspace : std::false_type{};

template<class T>
struct PossiblyAffineTrialColumnSubspace<
  T,
  mpl::enable_if_t<
    ::pressio::has_basis_matrix_typedef<T>::value
    && ::pressio::has_reduced_state_typedef<T>::value
    && ::pressio::has_full_state_typedef<T>::value
    /*
      requirements on the nested aliases
    */
    && pressio::is_vector_eigen<typename T::reduced_state_type>::value
    && std::is_copy_constructible< typename T::full_state_type>::value
    && std::is_copy_constructible< typename T::basis_matrix_type>::value
    && all_have_traits_and_same_scalar<
      typename T::reduced_state_type,
      typename T::full_state_type,
      typename T::basis_matrix_type>::value
    && Traits<typename T::full_state_type>::rank == 1
    && Traits<typename T::basis_matrix_type>::rank == 2
    /**/
    && std::is_copy_constructible<T>::value
    && !std::is_assignable<T&, T>::value
    && !std::is_assignable<T&, T&>::value
    /*
      requirements on the "create*" methods
    */
    && has_const_create_reduced_state_return_result<T>::value
    && has_const_create_full_state_return_result<T>::value
    && has_const_create_full_state_from_reduced_state<T>::value
    /*
      other
    */
    && has_const_map_from_reduced_state_return_void<T>::value
    && std::is_same<
      decltype( std::declval<T const>().basisOfTranslatedSpace() ),
      typename T::basis_matrix_type const &
      >::value
    && std::is_same<
      decltype(std::declval<T const>().translationVector()),
      const typename T::full_state_type &
      >::value
    && std::is_same<
        decltype( std::declval<T const>().basis() ),
        const typename T::basis_matrix_type &
      >::value
    && std::is_integral<
      decltype( std::declval<T const>().dimension() )
      >::value
    //
    && std::is_same<
      decltype( std::declval<T const>().isColumnSpace() ),
      bool
      >::value
    && std::is_same<
      decltype( std::declval<T const>().isRowSpace() ),
      bool
      >::value
    >
  > : std::true_type{};

}} //end namespace pressio::rom
#endif

#endif  // ROM_CONCEPTS_POSSIBLY_AFFINE_TRIAL_COLUMN_SUBSPACE_HPP_
