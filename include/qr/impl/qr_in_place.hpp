/*
//@HEADER
// ************************************************************************
//
// qr_in_place.hpp
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

#ifndef QR_IMPL_QR_IN_PLACE_HPP_
#define QR_IMPL_QR_IN_PLACE_HPP_

namespace pressio{ namespace qr{ namespace impl{

/* specialize for R_type == void */
template<
  typename matrix_type, typename algo,
  template <typename...> class Q_type
  >
class QRSolver<
  matrix_type, algo, true, void, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper<matrix_type>::value or
    containers::predicates::is_dense_matrix_wrapper<matrix_type>::value
    >
  > : public details::traits< QRSolver<matrix_type, algo,
				       true, void, Q_type>>::base_compute_t,
      public details::traits< QRSolver<matrix_type, algo,
				       true, void, Q_type>>::base_solve_t
{

  using this_t		= QRSolver<matrix_type, algo, true, void, Q_type>;
  using traits_t	= details::traits<this_t>;
  using base_compute_t	= typename traits_t::base_compute_t;
  using base_solve_t	= typename traits_t::base_solve_t;

  using impl_t	       = typename traits_t::impl_t;
  impl_t myImpl_;

  void computeThinImpl(matrix_type & A){
    myImpl_.computeThinInPlace(A);
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y) const {
    myImpl_.template doLinSolve<vector_t>(rhs, y);
  }

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  friend base_compute_t;
  friend base_solve_t;
};



// /* specialize for R_type != void */
// template<
//   typename matrix_type, typename algo,
//   typename R_type, template <typename...> class Q_type
//   >
// class QRSolver<
//   matrix_type, algo, true, R_type, Q_type,
//   ::pressio::mpl::enable_if_t<
//     meta::is_legitimate_r_type<R_type>::value and
//     containers::predicates::is_multi_vector_wrapper<matrix_type>::value
//     >
//   > : public details::traits< QRSolver<matrix_type, algo,
// 				       true, R_type, Q_type>>::base_compute_t,
//       public details::traits< QRSolver<matrix_type, algo,
// 				       true, R_type, Q_type>>::base_solve_t,
//       public details::traits< QRSolver<matrix_type, algo,
// 				       true, R_type, Q_type>>::base_Rfactor_t
// {

//   using this_t	       = QRSolver<matrix_type, algo, true, R_type, Q_type>;
//   using traits_t       = details::traits<this_t>;
//   using base_compute_t = typename traits_t::base_compute_t;
//   using base_solve_t   = typename traits_t::base_solve_t;
//   using base_Rfactor_t = typename traits_t::base_Rfactor_t;

//   using impl_t	       = typename traits_t::impl_t;
//   impl_t myImpl_       = {};

//   void computeThinImpl(matrix_type & A){
//     myImpl_.computeThinInPlace(A);
//   }

//   // here we have the Rfactor, maybe we should call the solver
//   // native to the Rfactor??? think
//   template <typename vector_t>
//   void solveImpl(const vector_t & rhs, vector_t & y) const{
//     myImpl_.template doLinSolve<vector_t, n>(rhs, y);
//   }

//   const R_type & cRefRFactorImpl() const {
//     return myImpl_.RFactor();
//   }

// public:
//   QRSolver() = default;
//   ~QRSolver() = default;

// private:
//   friend base_compute_t;
//   friend base_Rfactor_t;
//   friend base_solve_t;
// };


}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_QR_IN_PLACE_HPP_
