/*
//@HEADER
// ************************************************************************
//
// pressio_qr.hpp
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

#ifndef PRESSIO_QR_HPP_
#define PRESSIO_QR_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"

#include "qr/src/qr_ConfigDefs.hpp"
#include "qr/src/qr_fwd.hpp"
#include "qr/src/qr_algorithms_tags.hpp"

#include "qr/src/qr_meta.hpp"

#include "qr/src/base/qr_in_place_base.hpp"
#include "qr/src/base/qr_out_of_place_base.hpp"
#include "qr/src/base/qr_r_factor_base.hpp"
#include "qr/src/base/qr_solve_base.hpp"

#include "qr/src/qr_traits.hpp"

#include "qr/src/impl/qr_rfactor_solve_impl.hpp"
#include "qr/src/impl/qr_out_of_place.hpp"
#include "qr/src/impl/qr_in_place.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "qr/src/impl/eigen/qr_eigen_dense_out_of_place_impl.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "qr/src/impl/epetra/qr_epetra_multi_vector_tsqr_impl.hpp"
#include "qr/src/impl/epetra/qr_epetra_mv_householder_using_eigen_impl.hpp"
#include "qr/src/impl/epetra/qr_epetra_multi_vector_modified_gram_schmidt_impl.hpp"
#include "qr/src/impl/tpetra/qr_tpetra_multi_vector_tsqr_impl.hpp"
#include "qr/src/impl/tpetra/qr_tpetra_multi_vector_modified_gram_schmidt_impl.hpp"
#include "qr/src/impl/tpetra/qr_tpetra_mv_householder_using_eigen_impl.hpp"
#include "qr/src/impl/tpetra/qr_tpetra_block_multi_vector_tsqr_impl.hpp"
#endif

#endif
