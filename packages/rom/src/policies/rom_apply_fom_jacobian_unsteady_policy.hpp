/*
//@HEADER
// ************************************************************************
//
// rom_apply_fom_jacobian_unsteady_policy.hpp
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

#ifndef ROM_APPLY_FOM_JACOBIAN_UNSTEADY_HPP_
#define ROM_APPLY_FOM_JACOBIAN_UNSTEADY_HPP_

namespace pressio{ namespace rom{ namespace policy{

template <>
struct ApplyFomJacobianDefault<false>{

  //------------------------------------------
  // enabled for native c++
  //------------------------------------------
  template <
    typename fom_t, typename state_t,
    typename operand_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::pressio::containers::meta::is_array_pybind11<state_t>::value and
      !::pressio::containers::meta::is_array_pybind11<operand_t>::value
      > * = nullptr
#endif
    >
  auto evaluate(const fom_t	& fomObj,
		const state_t   & yFOM,
		const operand_t & B,
		time_t		  t) const
    -> decltype(fomObj.applyJacobian(*yFOM.data(), *B.data(), t))
  {
    return fomObj.applyJacobian(*yFOM.data(), *B.data(), t);
  }

  template <
    typename fom_t, typename state_t, typename operand_t,
    typename result_t, typename time_t
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      mpl::not_same<fom_t, pybind11::object>::value and
      !::pressio::containers::meta::is_array_pybind11<state_t>::value and
      !::pressio::containers::meta::is_array_pybind11<operand_t>::value
      > * = nullptr
#endif
    >
  void evaluate(const fom_t	  & fomObj,
		const state_t	  & yFOM,
		const operand_t & B,
		result_t	  & out,
		time_t		  t) const{
    fomObj.applyJacobian(*yFOM.data(), *B.data(), t, *out.data());
  }


#ifdef HAVE_PYBIND11
  //------------------------------------------
  // enabled when interfacing with python
  //------------------------------------------
  template <
    typename fom_t, typename state_t,
    typename operand_t, typename time_t
    , mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::pressio::containers::meta::is_array_pybind11<state_t>::value and
      ::pressio::containers::meta::is_array_pybind11<operand_t>::value and
      // we should have all data struct to be = pybind11::array_t
      mpl::is_same<state_t, operand_t>::value
      > * = nullptr
    >
  state_t evaluate(const fom_t	   & fomObj,
		   const state_t   & yFOM,
		   const operand_t & B,
		   time_t	     t) const{
    return fomObj.attr("applyJacobian")(yFOM, B, t);
  }

  template <
    typename fom_t, typename state_t, typename operand_t,
    typename result_t, typename time_t
    , mpl::enable_if_t<
      mpl::is_same<fom_t, pybind11::object>::value and
      ::pressio::containers::meta::is_array_pybind11<state_t>::value and
      ::pressio::containers::meta::is_array_pybind11<operand_t>::value and
      // because we should have all data struct to be = pybind11::array_t
      mpl::is_same<state_t, operand_t>::value
      > * = nullptr
    >
  void evaluate(const fom_t	  & fomObj,
		const state_t	  & yFOM,
		const operand_t & B,
		result_t	  & out,
		time_t		  t) const{
    fomObj.attr("applyJacobian")(yFOM, B, t, out);
  }
#endif

};

}}} //end namespace pressio::rom::policy
#endif
