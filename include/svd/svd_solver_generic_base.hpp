/*
//@HEADER
// ************************************************************************
//
// svd_solver_generic_base.hpp
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

#ifndef SVD_SVD_SOLVER_GENERIC_BASE_HPP_
#define SVD_SVD_SOLVER_GENERIC_BASE_HPP_


namespace pressio{ namespace svd{

template<typename derived_type>
class SolverBase
  : private utils::details::CrtpBase<SolverBase<derived_type>>
{

private:
  using matrix_t = typename svd::details::svd_traits<derived_type>::matrix_t;
  using sc_t = typename svd::details::svd_traits<derived_type>::scalar_t;
  using leftSvec_t = typename svd::details::svd_traits<derived_type>::lsv_t;
  using rightSvec_t = typename svd::details::svd_traits<derived_type>::rsv_t;
  using sval_t = typename svd::details::svd_traits<derived_type>::sval_t;

public:

  template<svdType svd_enum_value>
	::pressio::mpl::enable_if_t< svd_enum_value==svdType::truncated >
  compute(matrix_t & mat, int t, sc_t tol = 1e-12){
    this->underlying().template computeImpl<svd_enum_value>(mat, t, tol);
  }

  const leftSvec_t & cRefLeftSingularVectors() const {
    return this->underlying().cRefLeftSingularVectorsImpl();
  };

  const rightSvec_t & cRefRightSingularVectors() const {
    return this->underlying().cRefRightSingularVectorsImpl();
  };

  const sval_t & singularValues() const{
    return this->underlying().singularValuesImpl();
  };

private:
  SolverBase() = default;
  ~SolverBase() = default;

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<SolverBase<derived_type>>;

};//end class

} // end namespace svd
}//end namespace pressio
#endif  // SVD_SVD_SOLVER_GENERIC_BASE_HPP_
