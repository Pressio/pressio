/*
//@HEADER
// ************************************************************************
//
// pressio_containers.hpp
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

#ifndef PRESSIO_CONTAINERS_HPP_
#define PRESSIO_CONTAINERS_HPP_

/*************** IMPORTANT **************
The order below matters because this is how we make sure
headers are included (and classes found) in the proper order.
****************************************/

// dependencies
#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"

// generic
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#endif
#include "containers/src/containers_fwd.hpp"
#include "containers/src/containers_wrapped_types_enum.hpp"
#include "containers/src/containers_shared_traits.hpp"

// common predicates
#include "containers/src/predicates/typedefs/containers_has_communicator_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_data_map_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_global_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_local_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_scalar_typedef.hpp"
#include "containers/src/predicates/containers_has_method_size.hpp"
#include "containers/src/predicates/containers_has_method_extent.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers/src/predicates/containers_is_teuchos_rcp.hpp"
#endif

// include all detection for native types
#include "containers/src/predicates/native_types_detection/containers_native_all.hpp"

// concrete Vector/Matrix/MV/Tensor classes
#include "containers/src/vector/pressio_containers_vector_include.hpp"
#include "containers/src/dense_matrix/pressio_containers_dense_matrix_include.hpp"
#include "containers/src/sparse_matrix/pressio_containers_sparse_matrix_include.hpp"
#include "containers/src/multi_vector/pressio_containers_multi_vector_include.hpp"
#include "containers/src/tensor/pressio_containers_tensor_include.hpp"

// expressions
#include "containers/src/expressions/pressio_containers_expressions_include.hpp"

// experimental
#include "containers/src/experimental/containers_multivector_set.hpp"

// others
#include "containers/src/predicates/containers_is_wrapper.hpp"
#include "containers/src/predicates/containers_not_wrapper.hpp"
#include "containers/src/predicates/containers_are_wrappers.hpp"
#include "containers/src/predicates/containers_are_scalar_compatible.hpp"
#include "containers/src/collection/containers_static_collection.hpp"
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers/src/predicates/containers_have_matching_exe_space.hpp"
#endif

#endif
