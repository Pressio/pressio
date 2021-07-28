/*
//@HEADER
// ************************************************************************
//
// containers_dense_matrix_traits.hpp
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

#ifndef CONTAINERS_DENSE_MATRIX_CONTAINERS_DENSE_MATRIX_TRAITS_HPP_
#define CONTAINERS_DENSE_MATRIX_CONTAINERS_DENSE_MATRIX_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

/********************************
an arbitrary dense matrix
*******************************/
template <typename wrapped_type>
struct traits<
  DenseMatrix<wrapped_type>,
  mpl::enable_if_t<
    containers::predicates::is_admissible_as_dense_matrix_arbitrary<wrapped_type>::value
    >
  >
  : public containers_shared_traits<
  wrapped_type, WrappedPackageIdentifier::Arbitrary, false, 2>,
  public matrix_shared_traits<false>
{

  using scalar_t  = typename wrapped_type::value_type;
  using value_t   = typename wrapped_type::value_type;
  using size_t   = typename wrapped_type::size_type;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseArbitrary;
};


#ifdef PRESSIO_ENABLE_TPL_EIGEN
//***********************************
// eigen dense matrix
//***********************************
template <typename wrapped_type>
struct traits<
  DenseMatrix<wrapped_type>,
  mpl::enable_if_t<
    containers::predicates::is_dense_matrix_eigen<wrapped_type>::value
    >
  >
  : public containers_shared_traits<
  wrapped_type, WrappedPackageIdentifier::Eigen, true, 2>,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEigen;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr bool is_static = ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
                                      wrapped_type::ColsAtCompileTime != Eigen::Dynamic );
  static constexpr bool is_dynamic  = !is_static;

  using scalar_t  = typename wrapped_type::Scalar;
  using ordinal_t = typename wrapped_type::StorageIndex;
  using size_t    = ordinal_t;

  using subspan_ret_t = expressions::SubspanExpr<DenseMatrix<wrapped_type>>;
  using subspan_const_ret_t = expressions::SubspanExpr< const DenseMatrix<wrapped_type>>;
  using diag_ret_t = expressions::DiagExpr<DenseMatrix<wrapped_type>>;
  using diag_const_ret_t = expressions::DiagExpr< const DenseMatrix<wrapped_type>>;

};
#endif //PRESSIO_ENABLE_TPL_EIGEN

//**********************************
// for teuchos serial dense matrix
//**********************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename wrapped_type>
struct traits<
  DenseMatrix<wrapped_type>,
  mpl::enable_if_t<
    containers::predicates::is_dense_matrix_teuchos<wrapped_type>::value
    >
  >
  : public containers_shared_traits<
  wrapped_type, WrappedPackageIdentifier::Trilinos, true, 2
  >,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseTeuchosSerial;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_t = typename wrapped_type::scalarType;
  using ordinal_t = typename wrapped_type::ordinalType;
  using size_t    = ordinal_t;

  // for now, this must be empty until we enable support for subspanning a kokkos matrix
  using subspan_ret_t = void;
  using subspan_const_ret_t = void;
};
#endif

//***********************************
// epetra dense distributed matrix
//***********************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename wrapped_type>
struct traits<
  DenseMatrix<wrapped_type>,
  mpl::enable_if_t<
    containers::predicates::is_admissible_as_dense_matrix_epetra<wrapped_type>::value
    >
  >
  : public containers_shared_traits<
  wrapped_type, WrappedPackageIdentifier::Trilinos, false, 2>,
    public matrix_shared_traits<false>
{
  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEpetra;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_t = double;
  using local_ordinal_t = int;
  using global_ordinal_t = int;
  using size_t    = global_ordinal_t;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_BlockMap;
};
#endif

//*******************************
// Kokkos dense matrix
//*******************************
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename wrapped_type>
struct traits<
  DenseMatrix<wrapped_type>,
    ::pressio::mpl::enable_if_t<
      containers::predicates::is_dense_matrix_kokkos<wrapped_type>::value
    >
  >
  : public containers_shared_traits<
  wrapped_type, WrappedPackageIdentifier::Kokkos,
  true, //true because kokkos is for shared mem
  2
  >,
  public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
      wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseKokkos;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  using scalar_t	  = typename wrapped_type::traits::value_type;
  using layout		  = typename wrapped_type::traits::array_layout;
  using ordinal_t	  = typename wrapped_type::traits::size_type;
  using size_t		  = ordinal_t;

  using execution_space   = typename wrapped_type::traits::execution_space;
  using memory_space	  = typename wrapped_type::traits::memory_space;
  using device_t	  = typename wrapped_type::traits::device_type;
  using memory_traits	  = typename wrapped_type::traits::memory_traits;
  using host_mirror_space = typename wrapped_type::traits::host_mirror_space;
  using host_mirror_t     = typename wrapped_type::host_mirror_type;

  static constexpr bool is_static = wrapped_type::traits::rank_dynamic==0;
  static constexpr bool is_dynamic  = !is_static;

  using subspan_ret_t = expressions::SubspanExpr<DenseMatrix<wrapped_type>>;
  using subspan_const_ret_t = expressions::SubspanExpr< const DenseMatrix<wrapped_type>>;
  using diag_ret_t = expressions::DiagExpr<DenseMatrix<wrapped_type>>;
  using diag_const_ret_t = expressions::DiagExpr< const DenseMatrix<wrapped_type>>;

  static constexpr bool has_host_execution_space =
    (false
     #ifdef KOKKOS_ENABLE_SERIAL
     || std::is_same<execution_space, Kokkos::Serial>::value
     #endif
     #ifdef KOKKOS_ENABLE_OPENMP
     || std::is_same<execution_space, Kokkos::OpenMP>::value
     #endif
     );
};
#endif

// //*******************************
// // Pybind array
// //*******************************
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename wrapped_type>
// struct traits<
//   DenseMatrix<wrapped_type>,
//     mpl::enable_if_t<
//       containers::predicates::is_array_pybind<wrapped_type>::value
//     >
//   >
//   : public containers_shared_traits<
//   wrapped_type,
//   false, true, false,
//   WrappedPackageIdentifier::Pybind,
//   true, 2>,
//     public matrix_shared_traits<false>
// {

//   static constexpr WrappedMatrixIdentifier
//   wrapped_matrix_identifier = WrappedMatrixIdentifier::DensePybind;
//   static constexpr bool is_static = false;
//   static constexpr bool is_dynamic  = !is_static;

//   using scalar_t	 = typename wrapped_type::value_type;
//   using ordinal_t	 = std::size_t;
//   using size_t		 = ordinal_t;

//   using mut_proxy_t = decltype( std::declval<wrapped_type &>().mutable_unchecked() );
//   using proxy_t	    = decltype( std::declval<const wrapped_type &>().unchecked() );

//   using const_data_return_t = wrapped_type const *;
//   using data_return_t	    = wrapped_type *;
//   using reference_t	    = scalar_t &;
//   using const_reference_t   = scalar_t const &;
//   using diag_ret_t = expressions::DiagExpr<DenseMatrix<wrapped_type>>;
//   using diag_const_ret_t = expressions::DiagExpr< const DenseMatrix<wrapped_type>>;

//   // using span_ret_t	 = expressions::SpanExpr<Vector<wrapped_type>>;
//   // using span_const_ret_t = expressions::SpanExpr< const Vector<wrapped_type>>;
// };
// #endif

}}}//end namespace pressio::containers::details
#endif  // CONTAINERS_DENSE_MATRIX_CONTAINERS_DENSE_MATRIX_TRAITS_HPP_
