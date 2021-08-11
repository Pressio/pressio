/*
//@HEADER
// ************************************************************************
//
// rom_masked.hpp
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

#ifndef ROM_GALERKIN_IMPL_DECORATORS_ROM_MASKED_HPP_
#define ROM_GALERKIN_IMPL_DECORATORS_ROM_MASKED_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class masker_t, class maskable_t>
class Masked : public maskable_t
{
public:
  using maskable_data_type = typename maskable_t::data_type;
  using typename maskable_t::data_type;

private:
  std::reference_wrapper<const masker_t> masker_;
  mutable maskable_data_type maskedObj_;

public:
  Masked() = delete;
  Masked(const Masked &) = default;
  Masked & operator=(const Masked &) = default;
  Masked(Masked &&) = default;
  Masked & operator=(Masked &&) = default;
  ~Masked() = default;

  template <typename ... Args>
  Masked(const masker_t & maskerIn, Args && ... args)
    : maskable_t(std::forward<Args>(args)...),
      masker_(maskerIn),
      maskedObj_(maskerIn.createApplyMaskResult(*maskable_t::get().data()))
  {}

public:
  const maskable_data_type & get() const{
    return maskedObj_;
  }

  template<class galerkin_state_t, class fom_system_t, class scalar_t, typename ...Args>
  void compute(const galerkin_state_t & galerkinState,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & time,
	       Args && ...args) const
  {
    maskable_t::compute(galerkinState, fomSystemObj, time, std::forward<Args>(args)...);
    const auto & fullObj = maskable_t::get();
    masker_.get().applyMask(*fullObj.data(), time, *maskedObj_.data());
  }
};

}}}}
#endif  // ROM_GALERKIN_IMPL_DECORATORS_ROM_MASKED_HPP_
