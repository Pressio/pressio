
#ifndef ROMPP_APPS_LIN_ADV_DIFF_REACTION_2D_EPETRA_HPP_
#define ROMPP_APPS_LIN_ADV_DIFF_REACTION_2D_EPETRA_HPP_

#include "../../../CORE_ALL"
#include "Epetra_MpiComm.h"
#include <Epetra_config.h>
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include <cmath>
//#include "Tpetra_Experimental_BlockCrsMatrix_Helpers_def.hpp"
#include "Tpetra_Experimental_BlockCrsMatrix_Helpers_decl.hpp"


namespace rompp{ namespace apps{

class UnsteadyLinAdvDiffReac2dEpetra{
protected:
  using nativeVec	= Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix	= Epetra_CrsMatrix;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= nativeVec;
  using residual_type	= state_type;

public:
  UnsteadyLinAdvDiffReac2dEpetra(Epetra_MpiComm & comm,
				 int Nx, int Ny,
				 scalar_type K	  = 5.0,
				 scalar_type eps = 0.01)
    : comm_(comm), NxPhys_{Nx}, NyPhys_{Ny},
      Nx_{NxPhys_-2}, Ny_{NyPhys_}, // because Dirichlet in x, Neumann in y
      K_{K}, eps_{eps},
      dx_{Lx_/(Nx-1)},
      dy_{Ly_/(Ny-1)},
      dxSqInv_{1.0/(dx_*dx_)},
      dySqInv_{1.0/(dy_*dy_)},
      dx2Inv_{1./dx_},
      dy2Inv_{1./dy_}
  {}

private:
  void localIDToLiLj(int ID, int & li, int & lj) const{
    lj = ID/Nx_;
    li = ID % Nx_;
  }
  void globalIDToGiGj(int ID, int & gi, int & gj) const{
    gj = ID/Nx_;
    gi = ID % Nx_;
  }

public:
  state_type const & getInitialState() const{
    return *state_;
  };

  void createSolutionMap();
  void createGridMap();
  void setupPhysicalGrid();
  void setupFields();
  void setup();
  void assembleMatrix() const;
  void computeChem(const state_type &) const;

  int getNumGlobalDofs() const{ return numGlobalDof_; }
  int getNumLocalDofs() const{ return dofPerProc_; }

  std::shared_ptr<nativeVec>
  getX() const { return x_; }

  std::shared_ptr<nativeVec>
  getY() const { return y_; }

  std::shared_ptr<nativeMatrix>
  getMatrix() const { return A_; }

public:
  void residual(const state_type & yState,
		residual_type & rhs,
		scalar_type t) const{
    residual_impl(yState, rhs);
  }

  residual_type residual(const state_type & yState,
			 scalar_type t) const{
    residual_type R( *dofMap_ );
    residual_impl(yState, R);
    return R;
  };

protected:
  void residual_impl(const state_type & yState,
		     residual_type & R) const
  {
    this->assembleMatrix();
    this->computeChem(yState);
    A_->Multiply(false, yState, R);
    /* here we have R = A*state, where A = -conv+diffusion
     * so we need to sum the reaction part */
    R.Update(1., (*chemForcing_), 1.0);
  }

  void fillSource1(){
    for (auto i=0; i<gptPerProc_; i++){
	const scalar_type xij = (*x_)[i];
	const scalar_type yij = (*y_)[i];
	const scalar_type dx = (xij - oPtS1[0]);
	const scalar_type dy = (yij - oPtS1[1]);
	const scalar_type distance = dx*dx + dy*dy;
	(*s1_)[i] = ( std::sqrt(distance) <= rS1) ? 0.1 : 0.0;
    }
  }

  void fillSource2(){
    for (auto i=0; i<gptPerProc_; i++){
	const scalar_type xij = (*x_)[i];
	const scalar_type yij = (*y_)[i];
	const scalar_type dx = (xij - oPtS2[0]);
	const scalar_type dy = (yij - oPtS2[1]);
	const scalar_type distance = dx*dx + dy*dy;
	(*s2_)[i] = ( std::sqrt(distance) <= rS2) ? 0.1 : 0.0;
    }
  }

  void fillSource3(){
    s3_->PutScalar(0);
  }

protected:
  // radius where source 1 is active
  const scalar_type rS1 = {0.1};
  // radius where source 2 is active
  const scalar_type rS2 = {0.2};
  // center of the S1 source
  const std::array<scalar_type,2> oPtS1 = {{0.75, 1.2}};
  // center of the S2 source
  const std::array<scalar_type,2> oPtS2 = {{0.75, 1.}};

  const int numFields_ = 3;
  const int maxNonZeroPerRow_ = 7;
  scalar_type Lx_ = 1.0;
  scalar_type Ly_ = 2.0;
  std::array<scalar_type,2> xAxis_{0., 1.};
  std::array<scalar_type,2> yAxis_{0., 2.};

  Epetra_MpiComm & comm_;
  // physical grid points
  int NxPhys_{};
  int NyPhys_{};
  // we consider only inner grid along x, but
  // the full grid along y because of neumann BC
  // so actual grid points involved in calculations
  // is NOT same as physical grid
  int Nx_{};
  int Ny_{};

  scalar_type K_{};
  scalar_type eps_{};

  scalar_type dx_{};
  scalar_type dy_{};
  scalar_type dxSqInv_{};
  scalar_type dySqInv_{};
  scalar_type dx2Inv_{};
  scalar_type dy2Inv_{};

  // we need a map for the grid
  rcp<Epetra_Map> gridMap_;
  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = 3 * number of unknown grid points
  int numGlobalUnkGpt_;
  int *MyGlobalUnkGpt_;
  int gptPerProc_;

  // we need one map for the solution
  rcp<Epetra_Map> dofMap_;
  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = 3 * number of unknown grid points
  int numGlobalDof_;
  int *MyGlobalDof_;
  int dofPerProc_;

  rcp<nativeVec> x_;
  rcp<nativeVec> y_;
  rcp<nativeVec> u_;
  rcp<nativeVec> v_;
  mutable rcp<nativeMatrix> A_;
  mutable rcp<nativeVec> chemForcing_;
  mutable rcp<nativeVec> state_;

  mutable rcp<nativeVec> s1_;
  mutable rcp<nativeVec> s2_;
  mutable rcp<nativeVec> s3_;
};

}} //namespace rompp::apps
#endif
