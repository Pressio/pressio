/*
//@HEADER
// ************************************************************************
//
// apps_unsteady_nonlinear_adv_diff_reaction_flame_2d_sample_mesh_eigen.hpp
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

#ifndef PRESSIO_APPS_NONLIN_ADV_DIFF_REACTION_FLAME_2D_SAMPLE_MESH_EIGEN_HPP_
#define PRESSIO_APPS_NONLIN_ADV_DIFF_REACTION_FLAME_2D_SAMPLE_MESH_EIGEN_HPP_

#include "../../../CONTAINERS_ALL"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <cmath>
#include <array>

namespace pressio{ namespace apps{

class UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen{
protected:
  using this_t		= UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
  using nativeVec	= Eigen::VectorXd;

  using arr5_t		= std::array<int,5>;
  using arr2_t		= std::array<int,2>;

public:
  using scalar_type	= double;
  using state_type	= nativeVec;
  using velocity_type	= state_type;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using mv_t		= Eigen::MatrixXd;

  using graph_t		= std::vector<arr5_t>;
  using gids_map_t	= std::vector<arr2_t>;

  typedef Eigen::Triplet<scalar_type> Tr;
  using mat4_t		= Eigen::Matrix<scalar_type, 4, 4>;
  static constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();
  static constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  static constexpr auto two = ::pressio::utils::constants::two<scalar_type>();
  static constexpr auto three = ::pressio::utils::constants::three<scalar_type>();
  static constexpr auto four = ::pressio::utils::constants::four<scalar_type>();
  static constexpr auto oneHalf = one/two;
  static constexpr auto oneThird = one/three;
  static constexpr auto eight = two*four;
  static constexpr auto eightOvThree = eight/three;

public:
  UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen
  (int Nx, int Ny,
   const graph_t & graph,
   const gids_map_t & gidsMap,
   scalar_type K	= static_cast<scalar_type>(2), // cm^2/s
   scalar_type preExp	= static_cast<scalar_type>(5.5*1e12),
   scalar_type E	= static_cast<scalar_type>(8.0*1000.),
   std::array<scalar_type,3> W = {{2.016, 31.9, 18}})
    : K_{K},
      preExp_{preExp},
      negE_{-E},
      W_(W),
      rhoOvWH2{rho_/W_[0]},
      rhoOvWO2{rho_/W_[1]},
      WH2ovRho{W_[0]/rho_},
      WO2ovRho{W_[1]/rho_},
      WH2OovRho{W_[2]/rho_},
      Nx_{Nx},
      Ny_{Ny},
      graph_{graph},
      gidsMap_{gidsMap},
      dx_{Lx_/(Nx)},
      dy_{Ly_/(Ny)},
      dxSqInv_{one/(dx_*dx_)},
      dySqInv_{one/(dy_*dy_)},
      dx2Inv_{one/(two*dx_)},
      dy2Inv_{one/(two*dy_)}
  {}

public:
  state_type const & getInitialState() const{
    return state_;
  };

  void setupPhysicalGrid();
  void setupFields();
  void setup();

  // int getUnknownCount() const{ return this_t::numSpecies_*Nx_*Ny_; }
  scalar_type getDx() const { return dx_; };
  scalar_type getDy() const { return dy_; };
  nativeVec getX() const { return x_; }
  nativeVec getY() const { return y_; }
  // nativeVec getU() const { return u_; }
  // nativeVec getV() const { return v_; }

public:
  void velocity(const state_type & yState,
		scalar_type t, 
    velocity_type & rhs) const
  {
    velocity_impl(yState, rhs);
  }

  velocity_type velocity(const state_type & yState,
			 scalar_type t) const{
    // the velocity has size equal to the dofs of the cells where
    // we want velocity
    velocity_type R(numDof_r_);
    velocity_impl(yState, R);
    return R;
  };

  void jacobian(const state_type & u,  		
  		const scalar_type t,
      jacobian_type & jac) const{
    this->jacobian_impl(u, jac);
  }

  jacobian_type jacobian(const state_type & u,
  			 const scalar_type t) const{
    jacobian_type JJ(numDof_r_, numDof_);
    this->jacobian_impl(u, JJ);
    return JJ;
  }

  void applyJacobian(const state_type & yState,
  		     const mv_t & B,  		     
  		     scalar_type t, 
           mv_t & C) const{
    applyJacobian_impl(yState, B, C);
  }

  mv_t applyJacobian(const state_type & yState,
  		     const mv_t & B,
  		     scalar_type t) const{
    mv_t A( numDof_r_, B.cols() );
    A.setZero();
    applyJacobian_impl(yState, B, A);
    return A;
  };

private:
  void globalIDToGiGj(int ID, int & gi, int & gj) const{
    gj = ID/Nx_;
    gi = ID % Nx_;
  }

  void compute_dsdw(scalar_type wT,
		    scalar_type wH2,
		    scalar_type wO2,
		    scalar_type wH2O) const;

  void compute_sources(scalar_type wT,
		       scalar_type wH2,
		       scalar_type wO2,
		       scalar_type wH2O) const;

  void velocity_impl(const state_type & yState,
		     velocity_type & R) const;

  void jacobian_impl(const state_type & yState,
  		     jacobian_type & J) const;

  void applyJacobian_impl(const state_type & yState,
  			  const mv_t & B,
  			  mv_t & C) const{
    jacobian_type JJ(numDof_r_, numDof_);
    JJ.setZero();
    this->jacobian_impl(yState, JJ);
    C = JJ * B;
  }

protected:
  // T, H2, O2, H2O
  static constexpr int numSpecies_{4};
  const scalar_type Lx_{1.8};
  const scalar_type Ly_{0.9};
  const std::array<scalar_type,2> xAxis_{{0., Lx_}};
  const std::array<scalar_type,2> yAxis_{{0., Ly_}};

  const scalar_type initTemp_ = {300};
  const scalar_type rho_{0.00139};
  const scalar_type gasR_{8.314};
  const scalar_type Q_{9800};
  const scalar_type K_{};
  const scalar_type preExp_ = {};
  const scalar_type negE_ = {};
  const std::array<scalar_type,3> W_ = {};
  const scalar_type rhoOvWH2 = {};
  const scalar_type rhoOvWO2 = {};
  const scalar_type WH2ovRho = {};
  const scalar_type WO2ovRho = {};
  const scalar_type WH2OovRho = {};

  const std::array<scalar_type,numSpecies_> bcLeftGamma13_{{300, 0, 0, 0}};
  const std::array<scalar_type,numSpecies_> bcLeftGamma2_{{950, 0.0282, 0.2259, 0}};

  /*
    graph: contains a list such that
    1 0 3 2 -1

    first col: contains GIDs of cells where we want velocity
    1,2,3,4 col: contains GIDs of neighboring cells needed for stencil
		 the order of the neighbors is: east, north, west, south

    if a neighbor GID is = -1, it means that neighbor is outside of the domain boundary
    so if west GID is -1, it means that this velocity cell is near the left bundary
   */
  const graph_t & graph_;

  /*
    two columns list where:
    first col: GID of each cell enumerated according to the sample mesh
    second col: corresponding GID enumerated according to the full mesh
   */
  const gids_map_t & gidsMap_;

  // Nx_ and Ny_ contain the # of grid points for the FULL mesh
  const int Nx_{};
  const int Ny_{};
  const scalar_type dx_{};
  const scalar_type dy_{};
  const scalar_type dxSqInv_{};
  const scalar_type dySqInv_{};
  const scalar_type dx2Inv_{};
  const scalar_type dy2Inv_{};

  // auxiliary vars used in implementation
  mutable int gi_{};
  mutable int gj_{};

  mutable int cellGIDinFullMesh_{};
  mutable int cellGID_{};

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numSpecies_ * number_of_unknown_grid_points
  // _r_ stands for velocity
  int numGpt_;
  int numDof_;
  int numGpt_r_;
  int numDof_r_;

  nativeVec x_;
  nativeVec y_;
  nativeVec u_;
  nativeVec v_;
  mutable std::vector<Tr> tripletList;
  mutable jacobian_type J_;
  mutable nativeVec state_;

  // source term at each cell
  mutable std::array<scalar_type, numSpecies_> s_;

  // jacobian of source term: del s/del w
  mutable mat4_t dsdw_;

  // to label a grid point based on whether it belongs
  //to lower, mid or upper region of the flow, see diagram in paper
  mutable nativeVec regionLabel_;
};

}} //namespace pressio::apps
#endif
