/*
//@HEADER
// ************************************************************************
//
// rom_projected.hpp
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

#ifndef ROM_GALERKIN_IMPL_DECORATORS_ROM_PROJECTED_HPP_
#define ROM_GALERKIN_IMPL_DECORATORS_ROM_PROJECTED_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class ProjectorType, class SizeType, int Rank, class ProjectableType>
class Projected : public ProjectableType
{
public:
  Projected() = delete;
  Projected(const Projected &) = default;
  Projected & operator=(const Projected &) = default;
  Projected(Projected &&) = default;
  Projected & operator=(Projected &&) = default;
  ~Projected() = default;

  template <class ... Args, int _Rank = Rank, mpl::enable_if_t<_Rank==1,int> = 0>
  Projected(SizeType ext0,
	    const ProjectorType & projectorIn,
	    Args && ... args)
    : ProjectableType(std::forward<Args>(args)...),
      projector_(projectorIn),
      extents_{ext0}
  {}

  template <class ... Args, int _Rank = Rank, mpl::enable_if_t<_Rank==2,int> = 0>
  Projected(SizeType ext0, SizeType ext1,
	    const ProjectorType & projectorIn,
	    Args && ... args)
    : ProjectableType(std::forward<Args>(args)...),
      projector_(projectorIn),
      extents_{ext0, ext1}
  {}

public:
  template<class GalerkinOperatorType, int _Rank = Rank>
  mpl::enable_if_t<_Rank==1, GalerkinOperatorType>
  create() const
  {
    GalerkinOperatorType result;
    ::pressio::ops::resize(result, extents_[0]);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template<class GalerkinOperatorType, int _Rank = Rank>
  mpl::enable_if_t<_Rank==2, GalerkinOperatorType>
  create() const
  {
    GalerkinOperatorType result;
    ::pressio::ops::resize(result, extents_[0], extents_[1]);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template<class GalerkinOperatorType, class GalerkinStateType, class TimeType, class ...Args>
  void compute(GalerkinOperatorType & galerkinOperator,
	       const GalerkinStateType & galerkinState,
	       TimeType evalTime,
	       Args && ...args) const
  {
    ProjectableType::compute(galerkinState, evalTime, std::forward<Args>(args)...);
    // fomOperand = fom velocity or the fom apply jacobain object
    const auto & fomOperand = ProjectableType::get();

    // galerkinOperator is either galerkin rhs or jacobian
    projector_(fomOperand, evalTime, galerkinOperator);
  }

private:
  std::reference_wrapper<const ProjectorType> projector_;
  const std::array<SizeType, Rank> extents_ = {};
};

}}}}
#endif  // ROM_GALERKIN_IMPL_DECORATORS_ROM_PROJECTED_HPP_
