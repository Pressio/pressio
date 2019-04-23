
#ifndef ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_2D_BLOCK_TPETRA_HPP_
#define ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_2D_BLOCK_TPETRA_HPP_

#include "../../../CORE_ALL"
#include <Tpetra_Map.hpp>
#include <Tpetra_BlockVector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Details_DefaultTypes.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix_def.hpp>
#include <cmath>

namespace rompp{ namespace apps{

class UnsteadyNonLinAdvDiffReac2dBlockTpetra{
protected:
  template<typename T> using stdrcp = std::shared_ptr<T>;

  // use default types: associated with whatever
  // config was used to build trilinos/kokkos,
  // CUDA, OpenMP, threads, Serial.
  using map_t		= Tpetra::Map<>;
  using tpetVec		= Tpetra::Vector<>;
  using nativeVec	= Tpetra::Experimental::BlockVector<>;
  using nativeMV	= Tpetra::Experimental::BlockMultiVector<>;
  using nativeMatrix	= Tpetra::Experimental::BlockCrsMatrix<>;
  using graph_t		= typename nativeMatrix::crs_graph_type;

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
  using scalar_type	= ST;
  using state_type	= nativeVec;
  using residual_type	= state_type;

public:
  UnsteadyNonLinAdvDiffReac2dBlockTpetra(rcpcomm_t comm,
				      int Nx, int Ny,
				      scalar_type K   = 5.0,
				      scalar_type eps = 0.01)
    : K_{K}, eps_{eps},
      comm_(comm), rank_{comm_->getRank()},
      NxPhys_{Nx}, NyPhys_{Ny},
      // because Dirichlet in x, Neumann in y
      Nx_{NxPhys_-2}, Ny_{NyPhys_},
      dx_{Lx_/(Nx-1)},
      dy_{Ly_/(Ny-1)},
      dxSqInv_{1.0/(dx_*dx_)},
      dySqInv_{1.0/(dy_*dy_)},
      dx2Inv_{1./(2.*dx_)},
      dy2Inv_{1./(2.*dy_)}
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

  void assembleMatrix(const state_type & yState) const;
  void computeChem(const state_type &) const;

public:
  void residual(const state_type & yState,
		residual_type & rhs,
		scalar_type t) const{
    residual_impl(yState, rhs);
  }

  residual_type residual(const state_type & yState,
			 scalar_type t) const{
    residual_type R( *map_,
		     UnsteadyNonLinAdvDiffReac2dBlockTpetra::numSpecies_ );
    residual_impl(yState, R);
    return R;
  };

private:
  void localIDToLiLj(int ID, int & li, int & lj) const{
    lj = ID/Nx_;
    li = ID % Nx_;
  }
  void globalIDToGiGj(int ID, int & gi, int & gj) const{
    gj = ID/Nx_;
    gi = ID % Nx_;
  }

  void residual_impl(const state_type & yState,
		     residual_type & R) const
  {
    this->assembleMatrix(yState);
    this->computeChem(yState);
    A_->applyBlock(yState, R);
    /* here we have R = A*state, where A = -conv+diffusion
     * so we need to sum the reaction part */
    R.update(1., (*chemReac_), 1.0);
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
  const std::array<scalar_type,2> xAxis_{0., 1.};
  const std::array<scalar_type,2> yAxis_{0., 2.};

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

}} //namespace rompp::apps
#endif
