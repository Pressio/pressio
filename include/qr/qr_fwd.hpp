/*
//@HEADER
// ************************************************************************
//
// qr_fwd.hpp
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

#ifndef QR_QR_FWD_HPP_
#define QR_QR_FWD_HPP_

namespace pressio{  namespace qr{

// all fwd decl of mixin base classes
template<typename derived, typename matrix_t>
class QRInPlaceBase;

template<typename derived, typename matrix_t, typename Q_type>
class QROutOfPlaceBase;

template<typename derived, typename R_type>
class RFactorBase;

template<typename derived>
class QRSolveBase;


namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<
  typename matrix_type,
  typename R_t = void,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable= void>
class QRHouseholderDenseEigenMatrixWrapper;

template<
  typename matrix_type,
  typename R_t = void,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable= void>
class QRHouseholderEigenMultiVectorWrapper;
#endif


#if defined PRESSIO_ENABLE_TPL_TRILINOS
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<
  typename matrix_t,
  typename R_t,
  template <typename...> class Q_type = containers::MultiVector>
class EpetraMVHouseholderUsingEigen;

template<
  typename matrix_t,
  typename R_t,
  template <typename...> class Q_type = containers::MultiVector>
class TpetraMVHouseholderUsingEigen;
#endif

template<
  typename matrix_t,
  typename R_t,
  typename wrap_Q_type,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable = void>
class EpetraMVTSQR;

template<
  typename matrix_t,
  typename R_t,
  typename wrap_Q_type,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable = void>
class TpetraMVTSQR;

template<
  typename matrix_t,
  typename R_t,
  typename wrap_Q_type,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable = void>
class TpetraBlockMVTSQR;

template<
  typename matrix_t,
  typename R_t,
  typename wrap_Q_type = void,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable = void>
class ModGramSchmidtMVEpetra;

template<
  typename matrix_t,
  typename R_t,
  typename wrap_Q_type = void,
  template <typename...> class Q_type = containers::MultiVector,
  typename enable = void>
class ModGramSchmidtMVTpetra;

#endif //PRESSIO_ENABLE_TPL_TRILINOS


template<
  typename matrix_type,
  typename algorithm,
  bool in_place,
  typename R_type,
  template <typename...> class Q_type,
  typename enable = void>
class QRSolver;

}//end namespace pressio::qr::impl


template<
  typename matrix_type,
  typename algorithm,
  bool in_place = false,
  template <typename...> class Q_type = ::pressio::containers::MultiVector,
  typename enable = void
  >
using QRSolver = impl::QRSolver<matrix_type, algorithm, in_place, void, Q_type>;

template<
  typename matrix_type,
  typename algorithm,
  typename R_type,
  bool in_place = false,
  template <typename...> class Q_type = ::pressio::containers::MultiVector,
  typename enable = void
  >
using QRSolverWrapR = impl::QRSolver<matrix_type, algorithm, in_place, R_type, Q_type>;


}}//end namespace pressio::qr
#endif  // QR_QR_FWD_HPP_
