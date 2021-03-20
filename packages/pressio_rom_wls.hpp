/*
//@HEADER
// ************************************************************************
//
// pressio_rom_wls.hpp
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

#ifndef PRESSIO_ROM_WLS_HPP_
#define PRESSIO_ROM_WLS_HPP_

/*
   this header includes everything needed for WLS.
   NOTE that the order below matters!
   Includes are ordered in a logical way and this
   allows us to avoid ending up with a tangled system.
*/

// need all of the dependent packages
#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"
#include "pressio_optimizers.hpp"
#include "pressio_solvers.hpp"
#include "pressio_ode.hpp"

// common classes for rom
#include "rom/src/pressio_rom_common.hpp"

// wls classes
#include "rom/src/utils/rom_utils_set_gen_coordinates.hpp"
#include "rom/src/wls/rom_wls_types.hpp"
#include "rom/src/wls/rom_wls_jacobian_updating_tag.hpp"
#include "rom/src/wls/rom_wls_jacobians_container.hpp"
#include "rom/src/wls/rom_wls_preconditioners.hpp"

#include "rom/src/wls/predicates/rom_wls_is_legitimate_preconditioner_type.hpp"
#include "rom/src/wls/predicates/rom_wls_is_legitimate_jacobian_updating_tag.hpp"
#include "rom/src/wls/time_schemes/rom_wls_implicit_euler.hpp"
#include "rom/src/wls/time_schemes/rom_wls_bdf2.hpp"
#include "rom/src/wls/time_schemes/rom_wls_select_timescheme_helper.hpp"

#include "rom/src/wls/rom_wls_hessian_gradient_system_api.hpp"
#include "rom/src/wls/rom_wls_hessian_and_gradient_sequential_policy.hpp"
#include "rom/src/wls/rom_wls_solve_windows.hpp"

#endif
