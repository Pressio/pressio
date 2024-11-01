/*
//@HEADER
// ************************************************************************
//
// solvers_linear_eigen_iterative_matrix_free_impl.hpp
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

#ifndef PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_EIGEN_ITERATIVE_MATRIX_FREE_IMPL_HPP_
#define PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_EIGEN_ITERATIVE_MATRIX_FREE_IMPL_HPP_

#include "solvers_linear_iterative_base.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace pressio { namespace linearsolvers{

template<typename UserDefinedOperatorType>
class OperatorWrapper;

}}

namespace Eigen {
  namespace internal {

    template<typename UserDefinedOperatorType>
    struct traits< pressio::linearsolvers::OperatorWrapper<UserDefinedOperatorType> >
      : public Eigen::internal::traits<
      Eigen::Matrix<typename UserDefinedOperatorType::scalar_type,-1,-1>
      >
    {};

  }
}

namespace pressio { namespace linearsolvers{

template<typename UserDefinedOperatorType>
class OperatorWrapper :
    public Eigen::EigenBase<
  OperatorWrapper<UserDefinedOperatorType>
  >
{
public:
  using Scalar = typename UserDefinedOperatorType::scalar_type;
  using RealScalar = Scalar;
  using StorageIndex = int;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
  };

  OperatorWrapper() = default;

  OperatorWrapper(UserDefinedOperatorType const & valueIn)
    : m_userOperator(&valueIn) {}


  int rows() const { return m_userOperator->rows(); }
  int cols() const { return m_userOperator->cols(); }

  template<typename Rhs>
  Eigen::Product<OperatorWrapper<UserDefinedOperatorType>, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs>& x) const{
    using r_t = Eigen::Product<OperatorWrapper<UserDefinedOperatorType>, Rhs, Eigen::AliasFreeProduct>;
    return r_t(*this, x.derived());
  }

  void replace(const UserDefinedOperatorType & opIn) {
    m_userOperator = &opIn;
  }

  template<class OperandT, class ResultT>
  void applyAndAddTo(OperandT const & operand, ResultT & out) const {
    //
    // compute: out += operator * operand

    // out += *m_userOperator * operand;
    m_userOperator->applyAndAddTo(operand, out);
  }

private:
  UserDefinedOperatorType const *m_userOperator = nullptr;
};

}} // end namespace pressio::linearsolvers

namespace Eigen {
  namespace internal {

    template<typename Rhs, typename UserDefinedOpT>
    struct generic_product_impl<
      pressio::linearsolvers::OperatorWrapper<UserDefinedOpT>,
      Rhs, DenseShape, DenseShape, GemvProduct
      >
      : generic_product_impl_base<
         pressio::linearsolvers::OperatorWrapper<UserDefinedOpT>, Rhs,
         generic_product_impl<pressio::linearsolvers::OperatorWrapper<UserDefinedOpT>, Rhs>
      >
    {
      using Scalar = typename Product<
	pressio::linearsolvers::OperatorWrapper<UserDefinedOpT>,
	Rhs
	>::Scalar;

      template<typename Dest>
      static void scaleAndAddTo(
	 Dest& dst,
	 const pressio::linearsolvers::OperatorWrapper<UserDefinedOpT> & lhs,
	 const Rhs& rhs,
	 const Scalar& alpha)
      {
	// This method should implement "dst += alpha * lhs * rhs" inplace,
	// however, for iterative solvers, alpha is always equal to 1,
	// so let's not bother about it.
	assert(alpha==Scalar(1) && "scaling is not implemented");
	EIGEN_ONLY_USED_FOR_DEBUG(alpha);

	lhs.applyAndAddTo(rhs, dst);
      }
    };
  }
}

namespace pressio { namespace linearsolvers{ namespace impl{

template<typename TagType, typename UserDefinedOperatorType>
class EigenIterativeMatrixFree
  : public IterativeBase< EigenIterativeMatrixFree<TagType, UserDefinedOperatorType>>
{

public:
  using this_type = EigenIterative<TagType, UserDefinedOperatorType>;
  // using matrix_type	= UserDefinedOperatorType;
  using scalar_type = typename UserDefinedOperatorType::scalar_type;
  using solver_traits = ::pressio::linearsolvers::Traits<TagType>;
  using op_wrapper_t = OperatorWrapper<UserDefinedOperatorType>;
  using native_solver_type = typename solver_traits::template eigen_solver_type<op_wrapper_t>;
  using base_iterative_type = IterativeBase<this_type>;
  using iteration_type = typename base_iterative_type::iteration_type;

  static_assert(solver_traits::eigen_enabled == true,
   "the native solver must be from Eigen to use in EigenIterativeMatrixFree");
  static_assert(solver_traits::direct == false,
   "The native eigen solver must be iterative to use in EigenIterativeMatrixFree");

public:
  EigenIterativeMatrixFree() = default;

  iteration_type numIterationsExecuted() const{
    return mysolver_.iterations();
  }

  scalar_type finalError() const{
    return mysolver_.error();
  }

  void resetLinearSystem(const UserDefinedOperatorType& A)
  {
    mysolver_.setMaxIterations(this->maxIters_);
    m_wrapper.replace(A);
    mysolver_.compute(m_wrapper);
  }

  template <typename T>
  void solve(const T& b, T & y){
    mysolver_.setMaxIterations(this->maxIters_);
    y = mysolver_.solve(b);
  }

  template <typename T>
  void solve(const UserDefinedOperatorType & A, const T& b, T & y){
    this->resetLinearSystem(A);
    this->solve(b, y);
  }

private:
  friend base_iterative_type;
  native_solver_type mysolver_ = {};
  op_wrapper_t m_wrapper;
};

}}} // end namespace pressio::solvers::iterarive::impl
#endif
