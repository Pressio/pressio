/*
//@HEADER
// ************************************************************************
//
// apps_steady_linear_adv_diff_2d_epetra.hpp
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

#ifndef APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_HPP_
#define APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_HPP_

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

namespace pressio{ namespace apps{
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
//   using velocity_type	= state_type;

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
  {
    this->setup();
  }

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
  void createMap(){
    // total number of dof (we only consider the interior points)
    numGlobalDof_ = NxDof_ * NyDof_;
    myMap_ = std::make_shared<Epetra_Map>(numGlobalDof_, 0, comm_);
  }

  Epetra_Map const & getDataMap()const {
    return *myMap_;
  };

  void assembleMatrix() const
  {
    int gi{}; int gj{};
    int k{};
    std::array<int, 5> colInd{};
    std::array<scalar_type, 5> values{};
    scalar_type c_T_ip1{};
    scalar_type c_T_im1{};
    scalar_type c_T_ij {};
    scalar_type c_T_jp1{};
    scalar_type c_T_jm1{};
    int gid_ip1{};
    int gid_im1{};
    int gid_jp1{};
    int gid_jm1{};

    for (auto i=0; i<dofPerProc_; i++)
    {
      // global ID of current point
      auto GID = MyGlobalDof_[i];
      // find global i,j of this point
      globalIDToGiGj(GID, gi, gj);

      const scalar_type uij = (*u_)[i];
      const scalar_type vij = (*v_)[i];

      k = 0;
      c_T_ij    = - invPrRe_*dxSqInv_ - invPrRe_*dySqInv_;
      values[k] = 2.0*c_T_ij;
      colInd[k] = GID;
      k++;

      if (gi>=1){
        c_T_im1 = invPrRe_*dxSqInv_ + uij*dx2Inv_;
        gid_im1 = GID - 1;
        values[k] = c_T_im1;
        colInd[k] = gid_im1;
        k++;
      }

      if (gi<NxDof_-1){
        c_T_ip1 = invPrRe_*dxSqInv_ - uij*dx2Inv_;
        gid_ip1 = GID + 1;
        values[k] = c_T_ip1;
        colInd[k] = gid_ip1;
        k++;
      }

      if (gj>=1 and gj<NyDof_-1){
        c_T_jm1 = invPrRe_*dySqInv_ + vij*dy2Inv_;
        gid_jm1 = GID - NxDof_;
        values[k] = c_T_jm1;
        colInd[k] = gid_jm1;
        k++;

        c_T_jp1 = invPrRe_*dySqInv_ - vij*dy2Inv_;
        gid_jp1 = GID + NxDof_;
        values[k] = c_T_jp1;
        colInd[k] = gid_jp1;
        k++;
      }

      if (gj==0){
        c_T_jp1 = 2.0*invPrRe_*dySqInv_;
        gid_jp1 = GID + NxDof_;
        values[k] = c_T_jp1;
        colInd[k] = gid_jp1;
        k++;
      }

      if (gj==NyDof_-1){
        c_T_jm1 = 2.0*invPrRe_*dySqInv_;
        gid_jm1 = GID - NxDof_;
        values[k] = c_T_jm1;
        colInd[k] = gid_jm1;
        k++;
      }

      if (A_->Filled())
        A_->ReplaceGlobalValues(GID, k,
              values.data(), colInd.data());
      else
        A_->InsertGlobalValues(GID, k,
             values.data(), colInd.data());
    }//loop

    if(!A_->Filled())
      A_->FillComplete();
  }

  void fillRhs() const{
    f_->PutScalar( static_cast<scalar_type>(0) );

    int gi{}; int gj{};
    for (auto i=0; i<dofPerProc_; i++)
    {
      // global ID of current point
      auto GID = MyGlobalDof_[i];
      // find global i,j of this point
      globalIDToGiGj(GID, gi, gj);

      if (gi==0){
        auto value = invPrRe_*dxSqInv_ + (*u_)[i]*dx2Inv_;
        (*f_)[i] = -value * Tleft_;
      }

      if (gi==NxDof_-1){
        auto value = invPrRe_*dxSqInv_ - (*u_)[i]*dx2Inv_;
        (*f_)[i] = -value * Tright_;
      }
    }//loop
  }

  void solve(){
    Epetra_LinearProblem Problem(A_.get(), T_.get(), f_.get());
    AztecOO Solver(Problem);
    Solver.Iterate(500, solveTolerance_);
    Solver.NumIters();
    Solver.TrueResidual();
  }

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
  void setup(){
    createMap();
    dofPerProc_= myMap_->NumMyElements();
    MyGlobalDof_ = myMap_->MyGlobalElements();

    x_ = std::make_shared<nativeVec>(*myMap_);
    y_ = std::make_shared<nativeVec>(*myMap_);
    u_ = std::make_shared<nativeVec>(*myMap_);
    v_ = std::make_shared<nativeVec>(*myMap_);
    int gi{}; int gj{};
    for (auto i = 0; i<dofPerProc_; i++){
      auto GID = MyGlobalDof_[i];
      globalIDToGiGj(GID, gi, gj);
      (*x_)[i] = dx_ + gi * dx_;
      (*y_)[i] = gj * dy_;
      auto xval = (*x_)[i];
      auto yval = (*y_)[i];
      (*u_)[i] = -std::sin(M_PI*xval) * std::cos(M_PI*yval);
      (*v_)[i] = std::cos(M_PI*xval) * std::sin(M_PI*yval);
    }

    A_ = std::make_shared<nativeMatrix>(Copy, *myMap_, maxNonZeroPerRow_);
    T_ = std::make_shared<nativeVec>(*myMap_);
    f_ = std::make_shared<nativeVec>(*myMap_);
  }

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

}} //namespace pressio::apps
#endif  // APPS_STEADY_LINEAR_ADV_DIFF2D_APPS_STEADY_LINEAR_ADV_DIFF_2D_EPETRA_HPP_
