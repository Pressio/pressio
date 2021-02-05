/*
//@HEADER
// ************************************************************************
//
// pressio_apps.hpp
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

#ifndef PRESSIO_APPS_HPP_
#define PRESSIO_APPS_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"

// burgers1d eigen
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "apps/src/burgers1d/apps_burgers1d_eigen.hpp"
#include "apps/src/burgers1d/apps_burgers1d_eigen_discrete_time_api.hpp"
#include "apps/src/swe2d/apps_swe2d_eigen.hpp"
#include "apps/src/swe2d/apps_swe2d_hyper_eigen.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "apps/src/burgers1d/apps_burgers1d_epetra.hpp"
#include "apps/src/burgers1d/apps_burgers1d_epetra_preconditioned.hpp"
#include "apps/src/burgers1d/apps_burgers1d_epetra_reduced_no_mask.hpp"
#include "apps/src/burgers1d/apps_burgers1d_tpetra.hpp"
#include "apps/src/burgers1d/apps_burgers1d_tpetra_block.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "apps/src/burgers1d/apps_burgers1d_kokkos.hpp"
#endif

#include "apps/src/burgers1d/arbitraryDataStructures/apps_burgers1d_arb_ds.hpp"
#include "apps/src/burgers1d/arbitraryDataStructures/apps_burgers1d_arb_ds_custom_dense_matrix.hpp"
#include "apps/src/burgers1d/arbitraryDataStructures/apps_burgers1d_arb_ds_custom_vector.hpp"
#include "apps/src/burgers1d/arbitraryDataStructures/apps_burgers1d_arb_ds_discrete_time_api_adapter.hpp"
#include "apps/src/burgers1d/arbitraryDataStructures/apps_burgers1d_arb_ds_continuous_time_api_adapter.hpp"

#include "apps/src/burgers1d/apps_burgers1d_gold_states_explicit.hpp"
#include "apps/src/burgers1d/apps_burgers1d_gold_states_implicit.hpp"

// steady 2d adv-diff
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "apps/src/steady_linear_adv_diff2d/apps_steady_linear_adv_diff_2d_epetra.hpp"
#include "apps/src/steady_linear_adv_diff2d/apps_steady_linear_adv_diff_2d_epetra_rom_adapter.hpp"
#endif

#endif
