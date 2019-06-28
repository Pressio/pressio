
#ifndef ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_2D_EIGEN_HPP_
#define ROMPP_APPS_NONLIN_ADV_DIFF_REACTION_2D_EIGEN_HPP_

#include "../../../ALGEBRA_ALL"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <array>

namespace rompp{ namespace apps{

class UnsteadyNonLinAdvDiffReac2dEigen{
protected:
  using this_t		= UnsteadyNonLinAdvDiffReac2dEigen;
  using nativeVec	= Eigen::VectorXd;
  using mv_t		= Eigen::MatrixXd;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= nativeVec;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

  typedef Eigen::Triplet<scalar_type> Tr;
  static constexpr auto zero = ::rompp::utils::constants::zero<scalar_type>();
  static constexpr auto one = ::rompp::utils::constants::one<scalar_type>();
  static constexpr auto two = ::rompp::utils::constants::two<scalar_type>();

public:
  UnsteadyNonLinAdvDiffReac2dEigen
  (int Nx, int Ny,
   scalar_type K   = static_cast<scalar_type>(5),
   scalar_type eps = static_cast<scalar_type>(0.01))
    : NxPhys_{Nx}, NyPhys_{Ny},
      Nx_{NxPhys_-2}, Ny_{NyPhys_}, // because Dirichlet in x, Neumann in y
      K_{K}, eps_{eps},
      dx_{Lx_/(Nx-1)},
      dy_{Ly_/(Ny-1)},
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
  void computeJacobian(const state_type &) const;

  int getUnknownCount() const{ return this_t::numSpecies_*Nx_*Ny_; }
  scalar_type getDiffusivity() const { return eps_; };
  scalar_type getReactionRate() const { return K_; };
  scalar_type getDx() const { return dx_; };
  scalar_type getDy() const { return dy_; };
  nativeVec getX() const { return x_; }
  nativeVec getY() const { return y_; }
  nativeVec getU() const { return u_; }
  nativeVec getV() const { return v_; }
  nativeVec getS1() const { return s1_; }
  nativeVec getS2() const { return s2_; }
  nativeVec getS3() const { return s3_; }

public:
  void residual(const state_type & yState,
		residual_type & rhs,
		scalar_type t) const{
    residual_impl(yState, rhs);
  }

  residual_type residual(const state_type & yState,
			 scalar_type t) const{
    residual_type R(numDof_);
    residual_impl(yState, R);
    return R;
  };

  void jacobian(const state_type & u,
		jacobian_type & jac,
		const scalar_type /*t*/) const{
    this->jacobian_impl(u, jac);
  }

  jacobian_type jacobian(const state_type & u,
			 const scalar_type t) const{
    jacobian_type JJ(u.size(), u.size());
    this->jacobian_impl(u, JJ);
    return JJ;
  }

  // computes: C = Jac B where B is a multivector
  void applyJacobian(const state_type & yState,
  		     const mv_t & B,
  		     mv_t & C,
  		     scalar_type t) const{
    applyJacobian_impl(yState, B, C);
  }

  mv_t applyJacobian(const state_type & yState,
  		     const mv_t & B,
  		     scalar_type t) const{
    mv_t A( yState.size(), B.cols() );
    applyJacobian_impl(yState, B, A);
    return A;
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
  void fillSource1();
  void fillSource2();
  void fillSource3();

  void residual_impl(const state_type & yState,
		     residual_type & R) const;

  void jacobian_impl(const state_type & yState,
		     jacobian_type & J) const;

  void applyJacobian_impl(const state_type & yState,
  			  const mv_t & B,
  			  mv_t & C) const{
    jacobian_type JJ(yState.size(), yState.size());
    JJ.setZero();
    this->jacobian_impl(yState, JJ);
    C = JJ * B;
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

  static constexpr int numSpecies_{3};
  const scalar_type Lx_{1.0};
  const scalar_type Ly_{2.0};
  const std::array<scalar_type,2> xAxis_{{0., 1.}};
  const std::array<scalar_type,2> yAxis_{{0., 2.}};

  // physical grid points
  const int NxPhys_{};
  const int NyPhys_{};

  // we consider only inner grid along x, but
  // the full grid along y because of neumann BC
  // so actual grid points involved in calculations
  // is NOT same as physical grid
  const int Nx_{};
  const int Ny_{};

  const scalar_type K_{};
  const scalar_type eps_{};

  const scalar_type dx_{};
  const scalar_type dy_{};
  const scalar_type dxSqInv_{};
  const scalar_type dySqInv_{};
  const scalar_type dx2Inv_{};
  const scalar_type dy2Inv_{};

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = 3 * number of unknown grid points
  int numGpt_;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = 3 * number of unknown grid points
  int numDof_;

  nativeVec x_;
  nativeVec y_;
  nativeVec u_;
  nativeVec v_;
  mutable std::vector<Tr> tripletList;
  mutable jacobian_type J_;
  mutable nativeVec state_;

  mutable nativeVec s1_;
  mutable nativeVec s2_;
  mutable nativeVec s3_;
};

}} //namespace rompp::apps
#endif
