
#ifndef ROMPP_APPS_STEADY_LIN_ADV_DIFF_2D_EPETRA_HPP_
#define ROMPP_APPS_STEADY_LIN_ADV_DIFF_2D_EPETRA_HPP_

#include "../../../CORE_ALL"
#include "Epetra_MpiComm.h"
#include <Epetra_config.h>
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO_config.h"
#include "AztecOO.h"
#include <cmath>

namespace rompp{ namespace apps{
/*
  See this:
http://demonstrations.wolfram.com/SteadyStateTwoDimensionalConvectionDiffusionEquation/

  Solve this:
  (Pr*Re)^-1 * (d2T/dx2 + d2T/dy2) - u dT/dx -v dT/dy = 0

  Domain (0,1) times (0,2)
  Dirichlet BC on left wall and right wall
  Tleft = 0 , Tright = 1
  homogeneous Neumann on top and bottom

  u = - sin(pi x) cos(pi y)
  v = cos(pi x) sin(pi y)

  Use 2nd order FD, with following grid:

  o x x x x o
  o x x x x o
  o x x x x o
  o x x x x o

  The x denotes an unknown. o : known values
  Enumeration of unknowns is done left-right from bottom to top
  as follows:

    ...
  o x x x x o

    4 5 6 7
  o x x x x o

    0 1 2 3
  o x x x x o

  MPI decomp is done (simply) by splitting blocks of the domains
  vertically, and ranks go from 0 to n_ranks-1 starting from bottom.
 */

class SteadyLinAdvDiff2dEpetra{
protected:
  using nativeVec	= Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using nativeMatrix	= Epetra_CrsMatrix;
  using scalar_type	= double;

// public:
//   /* these types exposed because need to be detected */
//   using scalar_type	= double;
//   using state_type	= Epetra_Vector;
//   using residual_type	= state_type;

public:
  SteadyLinAdvDiff2dEpetra(Epetra_MpiComm & comm,
			int Nx, int Ny,
			scalar_type Pr = 2.,
			scalar_type Re = 30.,
			scalar_type Tleft = 0.0,
			scalar_type Tright = 1.0)
    : comm_(comm), Nx_{Nx}, Ny_{Ny},
      NxDof_{Nx-2}, NyDof_{Ny}, // because Dirichlet in x, Neumann in y
      Pr_{Pr}, Re_{Re}, invPrRe_{1./(Pr_*Re_)},
      Tleft_{Tleft}, Tright_{Tright},
      dx_{Lx_/(Nx-1)},
      dy_{Ly_/(Ny-1)},
      dxSqInv_{1.0/(dx_*dx_)},
      dySqInv_{1.0/(dy_*dy_)},
      dx2Inv_{1./dx_},
      dy2Inv_{1./dy_}
  {}

private:
  void localIDToLiLj(int ID, int & li, int & lj) const{
    lj = ID/NxDof_;
    li = ID % NxDof_;
  }
  void globalIDToGiGj(int ID, int & gi, int & gj) const{
    gj = ID/NxDof_;
    gi = ID % NxDof_;
  }


public:
  void createMap();
  Epetra_Map const & getDataMap()const { return *myMap_; };
  void setup();
  void assembleMatrix() const;
  void fillRhs() const;
  void solve();

  void printStateToFile(std::string fileName){
    std::ofstream file;
    file.open( fileName );
    for(auto i=0; i < dofPerProc_; i++){
      file << std::fixed << std::setprecision(14) <<
	(*x_)[i] << " " << (*y_)[i] << " " << (*T_)[i];
      file << std::endl;
    }
    file.close();
  }

  int getNumGlobalDofs() const{ return numGlobalDof_; }
  int getNumLocalDofs() const{ return dofPerProc_; }

  std::shared_ptr<nativeVec>
  getX() const { return x_; }

  std::shared_ptr<nativeVec>
  getY() const { return y_; }

  std::shared_ptr<nativeVec>
  getState() const { return T_; }

  std::shared_ptr<nativeMatrix>
  getMatrix() const { return A_; }

  std::shared_ptr<nativeVec>
  getForcing() const { return f_; }

protected:
  const int maxNonZeroPerRow_ = 5;
  const scalar_type solveTolerance_ = 1e-12;
  scalar_type Lx_ = 1.0;
  scalar_type Ly_ = 2.0;
  std::array<scalar_type,2> xAxis_{{0., 1.}};
  std::array<scalar_type,2> yAxis_{{0., 2.}};

  Epetra_MpiComm & comm_;
  // physical grid points
  int Nx_{};
  int Ny_{};
  // dof's which are not same as physical
  // we consider inner grid along x, but
  // the full grid along y because of neumann BC
  int NxDof_{};
  int NyDof_{};

  scalar_type Pr_{2.};
  scalar_type Re_{30.};
  scalar_type invPrRe_{};
  scalar_type Tleft_{0.0};
  scalar_type Tright_{1.0};

  scalar_type dx_{};
  scalar_type dy_{};
  scalar_type dxSqInv_{};
  scalar_type dySqInv_{};
  scalar_type dx2Inv_{};
  scalar_type dy2Inv_{};

  rcp<Epetra_Map> myMap_;
  int numGlobalDof_;
  int *MyGlobalDof_;
  int dofPerProc_;

  rcp<nativeVec> x_;
  rcp<nativeVec> y_;
  rcp<nativeVec> u_;
  rcp<nativeVec> v_;
  mutable rcp<nativeMatrix> A_;
  mutable rcp<nativeVec> T_;
  mutable rcp<nativeVec> f_;

};

}} //namespace rompp::apps
#endif
