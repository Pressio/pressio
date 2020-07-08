/*
//@HEADER
// ************************************************************************
//
// rom_query_fom_velocity_unsteady.hpp
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

#ifndef ROM_QUERY_FOM_VELOCITY_UNSTEADY_HPP_
#define ROM_QUERY_FOM_VELOCITY_UNSTEADY_HPP_

namespace pressio{ namespace rom{

template <typename fom_t, typename state_t, typename rhs_t, typename time_t>
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  mpl::enable_if_t<
  mpl::not_same<fom_t, pybind11::object>::value and
  !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
  !::pressio::containers::predicates::is_vector_wrapper_pybind<rhs_t>::value,
  int > = 0
#else
  void
#endif
queryFomVelocityUnsteady(const fom_t & fomObj,
       const state_t & yFOM,
       rhs_t & rhs,
       const time_t & t)
{
  fomObj.velocity(*yFOM.data(), t, *rhs.data());
}

//------------------------------------------
// enabled for native c++
//------------------------------------------
// template <
//   typename fom_t, typename state_t, typename time_t
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   , mpl::enable_if_t<
//       mpl::not_same<fom_t, pybind11::object>::value and
//       !::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value,
//       int > = 0
// #endif
//   >
// auto queryFomVelocityUnsteady(const fom_t & fomObj,
//             const state_t & yFOM,
//             const time_t & t)
//   -> decltype(fomObj.velocity(*yFOM.data(), t))
// {
//   return fomObj.velocity(*yFOM.data(), t);
// }

// //------------------------------------------
// // enabled when interfacing with python
// //------------------------------------------
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <
//   typename state_t, typename rhs_t, typename time_t,
//   mpl::enable_if_t<
//     ::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value and
//     ::pressio::containers::predicates::is_vector_wrapper_pybind<rhs_t>::value,
//     int > = 0
//   >
// static void queryFomVelocityUnsteady(const pybind11::object & fomObj,
// 				     const state_t & yFOM,
// 				     rhs_t & rhs,
// 				     const time_t & t)
// {
//  *rhs.data() = fomObj.attr("velocity")(*yFOM.data(), t);
// }

// template <typename state_t, typename time_t>
// static mpl::enable_if_t<
//   ::pressio::containers::predicates::is_vector_wrapper_pybind<state_t>::value,
//   typename ::pressio::containers::details::traits<state_t>::wrapped_t
//   >
// queryFomVelocityUnsteady(const pybind11::object & fomObj,
// 			 const state_t & yFOM,
// 			 const time_t & t)
// {
//  return fomObj.attr("velocity")(*yFOM.data(), t);
// }
// #endif

}} //end namespace pressio::rom
#endif
