/*
//@HEADER
// ************************************************************************
//
// containers_wrapped_types_enum.hpp
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

#ifndef CONTAINERS_CONTAINERS_WRAPPED_TYPES_ENUM_HPP_
#define CONTAINERS_CONTAINERS_WRAPPED_TYPES_ENUM_HPP_

namespace pressio{ namespace containers{ namespace details {

/*--------------------------------------------
Wrapped library name
--------------------------------------------*/
enum class WrappedPackageIdentifier{
   Undefined,
   Eigen,
   Kokkos,
   Trilinos,
   Pybind,
   Arbitrary
};

/*--------------------------------------------
Identifier for wrapped vectors

Within a given package, like trilinos, we can
have multiple types of vectors, e.g. epetra, tpetra.
Same can be true for other packages.
--------------------------------------------*/
enum class WrappedVectorIdentifier{
   Undefined,
   EigenRowStatic,
   EigenColStatic,
   EigenRowDynamic,
   EigenColDynamic,
   KokkosDynamic,
   KokkosStatic,
   Epetra,
   Tpetra,
   TpetraBlock,
   TeuchosSerialDense,
   Arbitrary
};

/*--------------------------------------------
Identifier for wrapped matrix
--------------------------------------------*/
enum class WrappedMatrixIdentifier{
   Undefined,
   DenseEigen,
   SparseEigen,
   DenseKokkos,
   DenseEpetra,
   DenseTeuchosSerial,
   DenseArbitrary
};

/*--------------------------------------------
Identifier for wrapped multivector
--------------------------------------------*/
enum class WrappedMultiVectorIdentifier{
   Undefined,
   Eigen,
   Kokkos,
   Epetra,
   Tpetra,
   TpetraBlock,
   Arbitrary
};

/*--------------------------------------------
Identifier for wrapped tensor
--------------------------------------------*/
enum class WrappedTensorIdentifier{
   Undefined,
   Pybind,
   Arbitrary
};

}}} // end namespace pressio::containers::details
#endif  // CONTAINERS_CONTAINERS_WRAPPED_TYPES_ENUM_HPP_
