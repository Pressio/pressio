/*
//@HEADER
// ************************************************************************
//
// apps_unsteady_nonlinear_adv_diff_reaction_2d_block_tpetra.hpp
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

#ifndef PRESSIO_APPS_NONLIN_ADV_DIFF_REACTION_2D_BLOCK_TPETRA_HPP_
#define PRESSIO_APPS_NONLIN_ADV_DIFF_REACTION_2D_BLOCK_TPETRA_HPP_

#include "../apps_ConfigDefs.hpp"

#ifdef HAVE_TRILINOS
#include "../../../CONTAINERS_ALL"
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Experimental_BlockVector.hpp>
#include <Tpetra_Experimental_BlockMultiVector.hpp>
#include <Tpetra_Details_DefaultTypes.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

namespace pressio{ namespace apps{

class UnsteadyNonLinAdvDiffReac2dBlockTpetra{
protected:
  template<typename T> using stdrcp = std::shared_ptr<T>;

  using this_t		= UnsteadyNonLinAdvDiffReac2dBlockTpetra;

  // use default types: associated with whatever
  // config was used to build trilinos/kokkos,
  // i.e. CUDA, OpenMP, threads, Serial.
  using map_t		= Tpetra::Map<>;
  using tpetVec		= Tpetra::Vector<>;
  using nativeVec	= Tpetra::Experimental::BlockVector<>;
  using nativeMV	= Tpetra::Experimental::BlockMultiVector<>;
  using nativeMatrix	= Tpetra::Experimental::BlockCrsMatrix<>;
  using graph_t		= typename nativeMatrix::crs_graph_type;

  // typedefs for communicators
  using tcomm		= Teuchos::Comm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  // get types (does not matter which one you use to detect types,
  // vector,mv,matrix should all have the same types
  // since they are all default)
  using ST	= typename nativeVec::scalar_type;
  using LO	= typename nativeVec::local_ordinal_type;
  using GO	= typename nativeVec::global_ordinal_type;
  using NT	= typename nativeVec::node_type;

  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // from device we can get the device execution and memory space
  using device_t	= typename nativeVec::device_type;
  using mem_space_t	= typename device_t::memory_space;
  using exec_space_t	= typename device_t::execution_space;

  // just check execution space is same from default types
  static_assert
  (std::is_same<
   typename Tpetra::Details::DefaultTypes::execution_space,
   exec_space_t>::value, "");

public:
  // public typedefs exposed to be detected by PRESSIO
  using scalar_type	= ST;
  using state_type	= nativeVec;
  using velocity_type	= state_type;

public:
  UnsteadyNonLinAdvDiffReac2dBlockTpetra
  (rcpcomm_t comm,
   GO Nx, GO Ny,
   scalar_type K   = static_cast<scalar_type>(5),
   scalar_type eps = static_cast<scalar_type>(0.01))
    : K_{K}, eps_{eps},
      comm_(comm), rank_{comm_->getRank()},
      NxPhys_{Nx}, NyPhys_{Ny},
      // because Dirichlet in x, Neumann in y
      Nx_{NxPhys_-2}, Ny_{NyPhys_},
      dx_{Lx_/(Nx-1)},
      dy_{Ly_/(Ny-1)},
      dxSqInv_{utils::constants::one<ST>()/(dx_*dx_)},
      dySqInv_{utils::constants::one<ST>()/(dy_*dy_)},
      dx2Inv_{utils::constants::one<ST>()/(utils::constants::two<ST>()*dx_)},
      dy2Inv_{utils::constants::one<ST>()/(utils::constants::two<ST>()*dy_)}
  {}

public:
  state_type const & getInitialState() const{
    return *state_;
  };

  void createMap();
  void initGridAndVel();
  void initFields();
  void createGraph();
  void computeSource();
  void setup();
  void assembleFDMatrix() const;
  void computeChem(const state_type &) const;
  void computeJacobian(const state_type & yState) const;
  GO getUnknownCount() const{ return this_t::numSpecies_*Nx_*Ny_; }
  int getNumSpecies() const{ return this_t::numSpecies_; }
  rcpmap_t getDataMap() const { return map_; }
  rcpmap_t getPointMap() const {
    map_t pmap = state_->getPointMap();
    return Teuchos::rcp(new map_t(pmap));
  }

public:
  void velocity(const state_type & yState,		
		scalar_type t, 
    velocity_type & rhs) const{
    velocity_impl(yState, rhs);
  }

  velocity_type velocity(const state_type & yState,
			 scalar_type t) const{
    velocity_type R( *map_,
		     UnsteadyNonLinAdvDiffReac2dBlockTpetra::numSpecies_ );
    velocity_impl(yState, R);
    return R;
  };

  // computes: C = Jac B where B is a multivector
  void applyJacobian(const state_type & yState,
  		   const nativeMV & B,  		     
		     scalar_type t,
         nativeMV & C) const{
    applyJacobian_impl(yState, B, C);
  }

  nativeMV applyJacobian(const state_type & yState,
			 const nativeMV & B,
			 scalar_type t) const{
    nativeMV C( *map_, B.getBlockSize(), B.getNumVectors() );
    applyJacobian_impl(yState, B, C);
    return C;
  };

private:
  void globalIDToGiGj(GO ID, GO & gi, GO & gj) const{
    gj = ID/Nx_;
    gi = ID % Nx_;
  }

  void velocity_impl(const state_type & yState,
		     velocity_type & R) const{
    static constexpr auto zero = ::pressio::utils::constants::zero<ST>();
    static constexpr auto one = ::pressio::utils::constants::one<ST>();

    R.putScalar(zero);
    this->assembleFDMatrix();
    this->computeChem(yState);
    A_->applyBlock(yState, R);
    /* here we have R = A*state, where A = -conv+diffusion
     * so we need to sum the reaction part */
    R.update(one, (*chemReac_), one);
  }

  void applyJacobian_impl(const state_type & yState,
			  const nativeMV & B,
			  nativeMV & C) const{
    static constexpr auto zero = ::pressio::utils::constants::zero<ST>();
    C.putScalar(zero);
    computeJacobian(yState);
    A_->applyBlock(B, C);
  }

private:
  // num of species
  static constexpr int numSpecies_ = 3;

  // radius where source 1 is active
  const scalar_type rS1 = {0.1};
  // radius where source 2 is active
  const scalar_type rS2 = {0.2};
  // center of the S1 source
  const std::array<scalar_type,2> oPtS1 = {{0.75, 1.2}};
  // center of the S2 source
  const std::array<scalar_type,2> oPtS2 = {{0.75, 1.}};

  // reaction coefficient
  scalar_type K_{};
  // diffusivity
  scalar_type eps_{};

  // domain size
  const scalar_type Lx_ = 1.0;
  const scalar_type Ly_ = 2.0;
  const std::array<scalar_type,2> xAxis_{{0., 1.}};
  const std::array<scalar_type,2> yAxis_{{0., 2.}};

  // communicator
  rcpcomm_t comm_{};
  int rank_{};

  // FD matrix has max 5 nonzeros per row
  static constexpr int maxNonZeroPerRow_ = 5;

  // physical grid points (the full vertex-centered grid)
  GO NxPhys_{};
  GO NyPhys_{};
  // the unknowns are only at inner grid along x, but
  // the full grid along y because of neumann BC
  // so actual grid points involved in calculations
  // is NOT same as physical grid
  GO Nx_{};
  GO Ny_{};

  scalar_type dx_{};
  scalar_type dy_{};
  scalar_type dxSqInv_{};
  scalar_type dySqInv_{};
  scalar_type dx2Inv_{};
  scalar_type dy2Inv_{};

  // the graph needed to create BlockCrsMatrix
  stdrcp<graph_t> graph_{};

  // the map
  rcpmap_t map_{};
  // note that the degress of freedom for this problem
  // are = 3 * number of unknown grid points
  GO numGlobalGpt_{};
  std::vector<GO> myGIDs_{};
  GO gptPerProc_{};

  stdrcp<tpetVec> x_{};
  stdrcp<tpetVec> y_{};
  stdrcp<tpetVec> u_{};
  stdrcp<tpetVec> v_{};

  mutable stdrcp<nativeMatrix> A_{};
  mutable stdrcp<nativeVec> chemReac_{};
  mutable stdrcp<nativeVec> state_{};
  mutable stdrcp<nativeVec> source_{};
};

}} //namespace pressio::apps
#endif
#endif
