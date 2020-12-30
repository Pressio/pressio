/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_kokkos.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_KOKKOS_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_KOKKOS_HPP_

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

public:
  explicit Burgers1dKokkos(params_t params,
			   ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell),
      x_d_{"xd", Ncell_},
      x_h_{"xd", Ncell_},
      U_d_{"ud", Ncell_},
      U_h_{"ud", Ncell_}
  {
    this->setup();
  }

  Burgers1dKokkos() = delete;
  ~Burgers1dKokkos() = default;

public:
  state_type const getInitialState() const {
    return U_d_;
  };

  mv_d createApplyJacobianResult(const mv_d & B) const
  {
    mv_d A("AA", Ncell_, B.extent(1) );
    return A;
  }

  velocity_type createVelocity() const{
    velocity_type RR("RR", Ncell_);
    return RR;
  }

  void velocity(const state_type & u,
		const scalar_type /* t */,
		velocity_type & rhs) const
  {
    using func_t = VelocityFunctor<k1dLl_d, state_type, velocity_type, sc_t>;
    func_t F(mu_[0], mu_[1], mu_[2], Ncell_, dxInv_, x_d_, u, rhs);
    Kokkos::parallel_for(Ncell_, F);
  }

  jacobian_type createJacobian() const
  {
    const int numRows = Ncell_;
    const int numCols = Ncell_;
    const int numEnt = (Ncell_-1)*2 + 1;

    // state_type_h u_h = Kokkos::create_mirror_view(u);
    typename jacobian_type::row_map_type::non_const_type ptr ("ptr", numRows+1);
    {
      typename jacobian_type::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
      ptr_h[0] = 0;
      for (int lclRow = 0; lclRow < numRows; ++lclRow) {
	auto nElPerRow = (lclRow==0) ? 1 : 2;
	ptr_h[lclRow+1] = ptr_h[lclRow] + nElPerRow; // 1 entry in row 0, 2 in all others
      }
      Kokkos::deep_copy (ptr, ptr_h);
    }

    typename jacobian_type::index_type::non_const_type ind ("ind", numEnt);
    {
      typename jacobian_type::index_type::HostMirror ind_h = Kokkos::create_mirror_view (ind);
      ind_h[0] = 0;
      int k = 0;
      for (int lclRow = 1; lclRow < numRows; ++lclRow) {
	ind_h[++k] = lclRow-1;
	ind_h[++k] = lclRow;
      }
      Kokkos::deep_copy (ind, ind_h);
    }

    typename jacobian_type::values_type val ("val", numEnt);
    {
      typename jacobian_type::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
 //      val_h[0] = -dxInv_*u_h(0);
 //      int k = 0;
 //      for (int lclRow = 1; lclRow < numRows; ++lclRow) {
	// val_h[++k] = dxInv_*u_h[lclRow-1];
	// val_h[++k] = -dxInv_*u_h[lclRow];
 //      }
      Kokkos::deep_copy (val, val_h);
    }

    jacobian_type JJ("JJ", numRows, numCols, numEnt, val, ptr, ind);
    // this->jacobian(u, t, JJ);
    return JJ;
  }

  void jacobian(const state_type & u,
    const scalar_type t,
    jacobian_type & jac) const
  {
    // here the Jacobian is passed in as an argument.
    // to recompute it, use parallel for

    // state_type_h u_h = Kokkos::create_mirror_view(u);
    // jacobian_type JJ2("JJ2", numRows, numCols, numEnt, val, ptr, ind);
    // jac = JJ2;
    //Kokkos::deep_copy(jac.values, JJ2.values);

    using exe_space = typename Burgers1dKokkos::execution_space;

    using policy_type = Kokkos::RangePolicy<exe_space, int>;

    using func_t = JacobianFunctor<state_type, jacobian_type, sc_t>;
    func_t F(Ncell_, dxInv_, u, jac);
    Kokkos::parallel_for( policy_type(0, Ncell_), F);
    //Kokkos::parallel_for(Ncell_, F);
  }

  void applyJacobian(const state_type & y,
         const mv_d & B,
         scalar_type t,
         mv_d & A) const
  {
    auto JJ = createJacobian();
    jacobian(y,t,JJ);
    constexpr auto zero = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto one = ::pressio::utils::constants<sc_t>::one();

    const char ct = 'N';
    KokkosSparse::spmv(&ct, one, JJ, B, zero, A);
  }

private:
  void setup()
  {
    dx_ = (xR_ - xL_)/static_cast<sc_t>(Ncell_);
    dxInv_ = 1.0/dx_;

    // grid
    for (ui_t i=0; i<Ncell_; ++i){
      x_h_(i) = dx_*i + dx_*0.5;
      U_h_(i) = 1.0;
    }
    Kokkos::deep_copy(x_d_, x_h_);
    Kokkos::deep_copy(U_d_, U_h_);
  }

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
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_KOKKOS_HPP_
