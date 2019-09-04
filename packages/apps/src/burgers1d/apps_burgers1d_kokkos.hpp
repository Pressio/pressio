/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_kokkos.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef PRESSIOAPPS_BURGERS1D_KOKKOS_HPP_
#define PRESSIOAPPS_BURGERS1D_KOKKOS_HPP_

#include "../apps_ConfigDefs.hpp"

// this has to be here because it is seen after we include configDefs
#ifdef HAVE_KOKKOS

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include "apps_burgers1d_kokkos_functors.hpp"
#include "KokkosSparse_spmv.hpp"
#include <array>

namespace pressio{ namespace apps{

struct Burgers1dKokkos{

  using sc_t	= double;
  using params_t = std::array<sc_t, 3>;
  using ui_t = unsigned int;

  using klr = Kokkos::LayoutRight;
  using kll = Kokkos::LayoutLeft;

  // using k1dLr_d = Kokkos::View<sc_t*, klr, execution_space>;
  // using k1dLr_h = k1dLr_d::host_mirror_type;
  // {
  //   // host and device have same layout
  //   using h_layout = typename k1dLr_h::traits::array_layout;
  //   using d_layout = typename k1dLr_d::traits::array_layout;
  //   static_assert( std::is_same<h_layout, d_layout>::value,
  // 		   "Layout for d and h (mirrorw) view is not the same");
  // }
  using execution_space = Kokkos::DefaultExecutionSpace;

  using k1dLl_d = Kokkos::View<sc_t*, kll, execution_space>;
  using k1dLl_h = k1dLl_d::host_mirror_type;
  // host and device should have same layout
  using h_layout1 = typename k1dLl_h::traits::array_layout;
  using d_layout1 = typename k1dLl_d::traits::array_layout;
  static_assert( std::is_same<h_layout1, d_layout1>::value,
		 "Layout for d and h (mirrorw) 1d view is not the same");

  using k2dLl_d = Kokkos::View<sc_t**, kll, execution_space>;
  using k2dLl_h = k2dLl_d::host_mirror_type;
  // host and device have same layout
  using h_layout2 = typename k2dLl_h::traits::array_layout;
  using d_layout2 = typename k2dLl_d::traits::array_layout;
  static_assert( std::is_same<h_layout2, d_layout2>::value,
		 "Layout for d and h (mirrorw) 2d view is not the same");

  using ord_t = int;
  using crs_mat = KokkosSparse::CrsMatrix<sc_t, ord_t, execution_space, void, int>;

  using state_type_d	= k1dLl_d;
  using state_type_h	= k1dLl_h;
  using velocity_type_d	= state_type_d;
  using velocity_type_h	= state_type_h;
  using mv_d = k2dLl_d;
  using mv_h = k2dLl_h;

  using scalar_type	= sc_t;
  using state_type	= state_type_d;
  using velocity_type	= state_type_d;
  using jacobian_type   = crs_mat;
  using mvec_t	= mv_d;

public:
  explicit Burgers1dKokkos(params_t params,
			   ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell),
      x_d_{"xd", Ncell_},
      x_h_{"xd", Ncell_},
      U_d_{"ud", Ncell_},
      U_h_{"ud", Ncell_}{}

  Burgers1dKokkos() = delete;
  ~Burgers1dKokkos() = default;

public:
  void setup();

  state_type const getInitialState() const {
    return U_d_;
  };

  void velocity(const state_type & u,
		const scalar_type /* t */, 
    velocity_type & rhs) const;

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const;

  void applyJacobian(const state_type & y,
		     const mvec_t & B,
		     scalar_type t, 
         mvec_t & A) const;

  mvec_t applyJacobian(const state_type & y,
		       const mvec_t & B,
		       scalar_type t) const;

  void jacobian(const state_type & u,
		const scalar_type t, 
    jacobian_type & jac) const;

  jacobian_type jacobian(const state_type & u,
			 const scalar_type t) const;

private:
  // parameters
  params_t mu_;

  //x coord of left side of domain
  const sc_t xL_ = 0.0;
  //x coord of right side of domain
  const sc_t xR_ = 100.0;

  // # of cells
  ui_t Ncell_;
  // cell size
  sc_t dx_;
  // 1/dx
  sc_t dxInv_;

  // mesh on device
  k1dLl_d x_d_;
  // mesh on host
  k1dLl_h x_h_;

  // state on device
  mutable state_type U_d_;
  // state on host
  mutable state_type_h U_h_;

};//end class

}} //namespace pressio::apps
#endif
#endif
