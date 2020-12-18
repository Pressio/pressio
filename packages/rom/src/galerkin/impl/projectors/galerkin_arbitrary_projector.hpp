/*
//@HEADER
// ************************************************************************
//
// galerkin_arbitrary_projector.hpp
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

#ifndef ROM_GALERKIN_IMPL_PROJECTORS_GALERKIN_ARBITRARY_PROJECTOR_HPP_
#define ROM_GALERKIN_IMPL_PROJECTORS_GALERKIN_ARBITRARY_PROJECTOR_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

/*
  arbitrary projector computes "P^T * operand":

   - operand = fom velocity (needed for explicit and implicit stepping)
     recall indeed that:
	dot(x_rom) = P^T f(phi*x_rom)

   - operand = phi (needed for jacobian of residual for implicit stepping)
     recall indeed that for implicit time stepping we have:
	R = dot(x_rom) - P^t f(phi*x_rom)
	so dR/dx_rom = ... - P^T df/dx phi, and call B = df/dx phi

   so we need to compute below "P^T f" and "P^T B".
*/

/*
   user-defined ops
*/
template <typename matrix_type, typename ud_ops_type>
struct ArbitraryProjector
{
  static_assert
  (::pressio::containers::predicates::is_dense_matrix_wrapper<matrix_type>::value or
   ::pressio::containers::predicates::is_multi_vector_wrapper<matrix_type>::value,
   "For arbitrary projector the matrix type must be a wrapper");

  using matrixtraits = ::pressio::containers::details::traits<matrix_type>;
  using matrix_native_type = typename matrixtraits::wrapped_t;

  ArbitraryProjector() = delete;
  ArbitraryProjector(const ArbitraryProjector &) = default;
  ArbitraryProjector & operator=(const ArbitraryProjector &) = delete;
  ArbitraryProjector(ArbitraryProjector &&) = default;
  ArbitraryProjector & operator=(ArbitraryProjector &&) = delete;
  ~ArbitraryProjector() = default;

  ArbitraryProjector(const matrix_type & matIn, const ud_ops_type & udOps)
    : matrix_{matIn}, udOps_(&udOps){}

  ArbitraryProjector(matrix_type && matIn, const ud_ops_type & udOps)
    : matrix_(std::move(matIn)), udOps_(udOps){}

  ArbitraryProjector(matrix_native_type && matIn, const ud_ops_type & udOps)
    : matrix_(std::move(matIn)), udOps_(udOps){}

  ArbitraryProjector(const matrix_native_type & matIn, const ud_ops_type & udOps)
    : matrix_(matIn), udOps_(udOps){}

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<::pressio::rom::galerkin::concepts::galerkin_rhs<result_t>::value>
  apply(const ::pressio::containers::Vector<operand_t> & operand,
	result_t & result) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<result_t>::scalar_t;
    using cnst = ::pressio::utils::constants<scalar_t>;
    udOps_.get().product(::pressio::transpose(), cnst::one(),
			 *matrix_.data(), *operand.data(),
			 cnst::zero(), result);
  }

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<
    ::pressio::rom::galerkin::concepts::galerkin_jacobian<result_t>::value and
    (::pressio::containers::predicates::is_dense_matrix_wrapper<operand_t>::value or
     ::pressio::containers::predicates::is_multi_vector_wrapper<operand_t>::value)
    >
  apply(const operand_t & operand, result_t & result) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<result_t>::scalar_t;
    using cnst = ::pressio::utils::constants<scalar_t>;
    udOps_.get().product(::pressio::transpose(), ::pressio::nontranspose(),
			 cnst::one(), *matrix_.data(), *operand.data(),
			 cnst::zero(), result);
  }

private:
  const matrix_type matrix_;
  std::reference_wrapper<const ud_ops_type> udOps_;
};


/*
   void ops, use pressio ops
*/
template <typename matrix_type>
struct ArbitraryProjector<matrix_type, void>
{
  static_assert
  (::pressio::containers::predicates::is_dense_matrix_wrapper<matrix_type>::value or
   ::pressio::containers::predicates::is_multi_vector_wrapper<matrix_type>::value,
   "For arbitrary projector the matrix type must be a wrapper");

  using matrixtraits = ::pressio::containers::details::traits<matrix_type>;
  using matrix_native_type = typename matrixtraits::wrapped_t;
  using scalar_t = typename matrixtraits::scalar_t;
  using cnst = ::pressio::utils::constants<scalar_t>;

  ArbitraryProjector() = delete;
  ArbitraryProjector(const ArbitraryProjector &) = default;
  ArbitraryProjector & operator=(const ArbitraryProjector &) = delete;
  ArbitraryProjector(ArbitraryProjector &&) = default;
  ArbitraryProjector & operator=(ArbitraryProjector &&) = delete;
  ~ArbitraryProjector() = default;

  explicit ArbitraryProjector(const matrix_type & matIn)
    : matrix_{matIn}{}

  explicit ArbitraryProjector(matrix_type && matIn)
    : matrix_(std::move(matIn)){}

  explicit ArbitraryProjector(matrix_native_type && matIn)
    : matrix_(std::move(matIn)){}

  explicit ArbitraryProjector(const matrix_native_type & matIn)
    : matrix_(matIn){}

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<::pressio::rom::galerkin::concepts::galerkin_rhs<result_t>::value>
  apply(const ::pressio::containers::Vector<operand_t> & operand,
	result_t & result) const
  {
    assert(matrix_.extent(0) == operand.extent(0));
    assert(matrix_.extent(1) == result.extent(0));

    ::pressio::ops::product(::pressio::transpose(), cnst::one(),
			    matrix_, operand, cnst::zero(), result);
  }

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<
    ::pressio::rom::galerkin::concepts::galerkin_jacobian<result_t>::value and
    (::pressio::containers::predicates::is_dense_matrix_wrapper<operand_t>::value or
     ::pressio::containers::predicates::is_multi_vector_wrapper<operand_t>::value)
    >
  apply(const operand_t & operand, result_t & result) const
  {
    assert(matrix_.extent(0) == operand.extent(0));
    assert(matrix_.extent(1) == result.extent(0));
    assert(operand.extent(1) == result.extent(1));

    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
     			    cnst::one(), matrix_, operand, cnst::zero(), result);
  }

private:
  const matrix_type matrix_;
};

}}}}//end  namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_PROJECTORS_GALERKIN_ARBITRARY_PROJECTOR_HPP_
