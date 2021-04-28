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
template <typename data_type, typename ud_ops_type>
struct ArbitraryProjector
{
  static_assert
  (::pressio::containers::predicates::is_wrapper<data_type>::value,
   "For arbitrary projector the data type must be a wrapper");

  using traits = ::pressio::containers::details::traits<data_type>;
  using matrix_native_type = typename traits::wrapped_t;

  ArbitraryProjector() = delete;
  ArbitraryProjector(const ArbitraryProjector &) = default;
  ArbitraryProjector & operator=(const ArbitraryProjector &) = delete;
  ArbitraryProjector(ArbitraryProjector &&) = default;
  ArbitraryProjector & operator=(ArbitraryProjector &&) = delete;
  ~ArbitraryProjector() = default;

  ArbitraryProjector(const data_type & dataIn, const ud_ops_type & udOps)
    : data_{dataIn}, udOps_(&udOps){}

  ArbitraryProjector(data_type && dataIn, const ud_ops_type & udOps)
    : data_(std::move(dataIn)), udOps_(udOps){}

  ArbitraryProjector(matrix_native_type && dataIn, const ud_ops_type & udOps)
    : data_(std::move(dataIn)), udOps_(udOps){}

  ArbitraryProjector(const matrix_native_type & dataIn, const ud_ops_type & udOps)
    : data_(dataIn), udOps_(udOps){}

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<
    ::pressio::containers::details::traits<operand_t>::rank == 1 and
    ::pressio::containers::details::traits<result_t>::rank == 1 and
    (::pressio::rom::galerkin::constraints::velocity<result_t>::value or
     ::pressio::rom::galerkin::constraints::residual<result_t>::value)
    >
  apply(const operand_t & operand, result_t & result) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<result_t>::scalar_t;
    using cnst = ::pressio::utils::constants<scalar_t>;
    udOps_.get().product(::pressio::transpose(), cnst::one(),
			 *data_.data(), *operand.data(),
			 cnst::zero(), result);
  }

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<
    ::pressio::containers::details::traits<operand_t>::rank ==2 and
    ::pressio::containers::details::traits<result_t>::rank == 1 and
    ::pressio::rom::galerkin::constraints::galerkin_jacobian<result_t>::value
    >
  apply(const operand_t & operand, result_t & result) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<result_t>::scalar_t;
    using cnst = ::pressio::utils::constants<scalar_t>;
    udOps_.get().product(::pressio::transpose(), ::pressio::nontranspose(),
			 cnst::one(), *data_.data(), *operand.data(),
			 cnst::zero(), result);
  }

private:
  const data_type data_;
  std::reference_wrapper<const ud_ops_type> udOps_;
};


/*
   void ops, use pressio ops
*/
template <typename data_type>
struct ArbitraryProjector<data_type, void>
{
  static_assert
  (::pressio::containers::predicates::is_wrapper<data_type>::value,
   "For arbitrary projector the data type must be a wrapper");

  using traits = ::pressio::containers::details::traits<data_type>;
  using native_type = typename traits::wrapped_t;
  using scalar_t = typename traits::scalar_t;
  using cnst = ::pressio::utils::constants<scalar_t>;

  ArbitraryProjector() = delete;
  ArbitraryProjector(const ArbitraryProjector &) = default;
  ArbitraryProjector & operator=(const ArbitraryProjector &) = delete;
  ArbitraryProjector(ArbitraryProjector &&) = default;
  ArbitraryProjector & operator=(ArbitraryProjector &&) = delete;
  ~ArbitraryProjector() = default;

  explicit ArbitraryProjector(const data_type & dataIn)
    : data_{dataIn}{}

  explicit ArbitraryProjector(data_type && dataIn)
    : data_(std::move(dataIn)){}

  explicit ArbitraryProjector(native_type && dataIn)
    : data_(std::move(dataIn)){}

  explicit ArbitraryProjector(const native_type & dataIn)
    : data_(dataIn){}

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<::pressio::containers::details::traits<result_t>::rank == 1>
  apply(const operand_t & operand, result_t & result) const
  {
    ::pressio::ops::product(::pressio::transpose(), cnst::one(),
			    data_, operand, cnst::zero(), result);
  }

  template<typename operand_t, typename result_t>
  mpl::enable_if_t<::pressio::containers::details::traits<result_t>::rank >= 2>
  apply(const operand_t & operand, result_t & result) const
  {
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    cnst::one(), data_, operand,
			    cnst::zero(), result);
  }

  // template<typename operand_t, typename result_t>
  // mpl::enable_if_t<
  //   ::pressio::containers::details::traits<operand_t>::rank ==2 and
  //   ::pressio::rom::galerkin::constraints::galerkin_jacobian<result_t>::value
  //   >
  // apply(const operand_t & operand, result_t & result) const
  // {
  //   // ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
  //   //  			    cnst::one(), data_, operand, cnst::zero(), result);
  // }

private:
  const data_type data_;
};

}}}}//end  namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_PROJECTORS_GALERKIN_ARBITRARY_PROJECTOR_HPP_
