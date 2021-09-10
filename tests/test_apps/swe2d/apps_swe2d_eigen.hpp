
#ifndef PRESSIO_APPS_SWE_HPP_
#define PRESSIO_APPS_SWE_HPP_

#include "Eigen/Core"
#include "rusanovflux.hpp"
#include <chrono>

namespace pressio{ namespace apps{

class swe2d
{
  using eigVec = Eigen::VectorXd;
  using ui_t = unsigned int;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using eig_sp_mat = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using jacobian_type	= eig_sp_mat;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  swe2d() {
    this->setup();
  }

  swe2d(const scalar_type Lx,
	const scalar_type Ly,
        const int nx,
	const int ny,
	const scalar_type params[3])
    : nx_(nx), ny_(ny),
      Lx_(Lx), Ly_(Ly),
      g_(params[0]),
      mu_ic_(params[1]),
      mu_f_(params[2])
  {
    this->setup();
  }

  swe2d(const scalar_type Lx,
	const scalar_type Ly,
        const int nx,
	const int ny)
    : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly)
  {
    this->setup();
  }

  swe2d(const int n) : nx_(n), ny_(n)
  {
    this->setup();
  }

  swe2d(const int n, const scalar_type params[3])
    : nx_(n), ny_(n),
      g_(params[0]),
      mu_ic_(params[1]),
      mu_f_(params[2])
  {
    this->setup();
  }

  int numDofs() const{
    return nDofs_;
  }

  state_type getGaussianIC(scalar_type mu) const
  {
    state_type U0(3*nx_*ny_);
    for (int j=0; j < ny_; j++)
    {
      for (int i=0; i < nx_; i++)
      {
        auto sid = state_index_mapper(0,i,j);
        U0(sid) = 1. + mu_ic_*exp( - ( pow( xGrid_(index_mapper(i,j)) - 1.5,2)
				       + pow( yGrid_(index_mapper(i,j)) - 1.5,2) ));
        U0(sid+1) = 0.;
        U0(sid+2) = 0.;
      }
    }
    return U0;
  }

  void setParams(const scalar_type params[3])
  {
    g_ = params[0];
    mu_ic_ = params[1];
    mu_f_ = params[2];
  }

  velocity_type createVelocity() const {
    velocity_type V(3*nx_*ny_);
    return V;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3*nx_*ny_,3*nx_*ny_);
    return J;
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
  {
    jacobian_type JA(3*nx_*ny_, A.cols() );
    return JA;
  }


  void velocity(const state_type & U,
                const scalar_type /*t*/,
                velocity_type & V) const
  {
    scalar_type FL[3];
    scalar_type FR[3];
    scalar_type FU[3];
    scalar_type FD[3];
    scalar_type forcing[3];
    const std::array<scalar_type, 2> nx = { 1, 0 };
    const std::array<scalar_type, 2> ny = { 0, 1 };

    int isL;
    int isR;
    int isU;
    int isD;
    for (int j=0; j < ny_ ; j++){
      for (int i=0; i < nx_  ; i++){
        isL = state_index_mapper(0,i-1,j);
        isR = state_index_mapper(0,i,j);
	rusanovFluxFullStateIn(FL,U,isL,isR,nx,g_);
        isL = state_index_mapper(0,i,j);
        isR = state_index_mapper(0,i+1,j);
	rusanovFluxFullStateIn(FR,U,isL,isR,nx,g_);

        isD = state_index_mapper(0,i,j-1);
        isU = state_index_mapper(0,i,j);
	rusanovFluxFullStateIn(FD,U,isD,isU,ny,g_);

        isD = state_index_mapper(0,i,j);
        isU = state_index_mapper(0,i,j+1);
	rusanovFluxFullStateIn(FU,U,isD,isU,ny,g_);

        auto is = state_index_mapper(0,i,j);
        forcing[0] = 0;
        forcing[1] = mu_f_*U(is+2);
        forcing[2] = mu_f_*U(is+1);
	for (int k=0;k<3;k++){
	  V(is+k) = -1./dx_*(FR[k] - FL[k]) - 1./dy_*(FU[k] - FD[k]) + forcing[k];
	}
      }
    }
  }

  void jacobian(const state_type & U,
                const scalar_type /*t*/,
                jacobian_type & jac) const
    {
      scalar_type JL_L[3][3];
      scalar_type JR_L[3][3];
      scalar_type JL_R[3][3];
      scalar_type JR_R[3][3];
      scalar_type JD_D[3][3];
      scalar_type JU_D[3][3];
      scalar_type JD_U[3][3];
      scalar_type JU_U[3][3];

      std::array<scalar_type, 2> nx = { 1, 0 };
      std::array<scalar_type, 2> ny = { 0, 1 };

      tripletList.clear();
      for (int sid=0; sid < nx_*ny_  ; sid++){
        const auto ij = get_ij_from_gid(sid);
        const auto gid = state_index_mapper(0,ij[0],ij[1]);
        const auto gid_im1 = state_index_mapper(0,ij[0] - 1,ij[1]);
        const auto gid_ip1 = state_index_mapper(0,ij[0] + 1,ij[1]);
        const auto gid_jm1 = state_index_mapper(0,ij[0],ij[1] - 1);
        const auto gid_jp1 = state_index_mapper(0,ij[0],ij[1] + 1);

        rusanovFluxJacobianFullStateIn(JL_L,JR_L,U,gid_im1,gid,nx,g_);
        rusanovFluxJacobianFullStateIn(JL_R,JR_R,U,gid,gid_ip1,nx,g_);
        rusanovFluxJacobianFullStateIn(JD_D,JU_D,U,gid_jm1,gid,ny,g_);
        rusanovFluxJacobianFullStateIn(JD_U,JU_U,U,gid,gid_jp1,ny,g_);
        for (int j = 0; j < 3; ++j) {
          for (int i = 0; i < 3; ++i) {
            auto val = -1./dx_*(JL_R[i][j] - JR_L[i][j]) - 1./dy_*(JD_U[i][j] - JU_D[i][j]);
            tripletList.push_back( Eigen::Triplet<scalar_type>( gid+i,gid+j,val ) );

            val =  1./dx_*JL_L[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_type>( gid+i,gid_im1+j,val ) );

            val = -1./dx_*JR_R[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_type>( gid+i,gid_ip1+j,val ) );

            val = 1./dy_*JD_D[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_type>( gid+i,gid_jm1+j,val ) );

            val = -1./dy_*JU_U[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_type>( gid+i,gid_jp1+j,val ) );

          }
        }
        // add forcing terms
        tripletList.push_back( Eigen::Triplet<scalar_type>( gid+1,gid+2,mu_f_ ) );
        tripletList.push_back( Eigen::Triplet<scalar_type>( gid+2,gid+1,mu_f_ ) );
      }
      jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void applyJacobian(const state_type &U,
		     const dense_matrix_type & A,
		     scalar_type t,
		     dense_matrix_type &JA) const
  {
    jacobian(U,t,jac_);
    JA = jac_*A;
  }

  // finite difference Jacobian, keep for validation
  void jacobian_fd(const state_type & U,
                   const scalar_type /*t*/,
                   jacobian_type & jac) const
  {
    velocity_type V0(nDofs_);
    velocity_type Vp(nDofs_);
    state_type Up(U);
    velocity(U,0.,V0);
    const scalar_type eps = 1e-5;
    tripletList.clear();
    for (int i = 0; i < nDofs_; i++){
      Up(i) += eps;
      velocity(Up,0.,Vp);
      auto colVal = (Vp - V0)/eps;
      Up(i) -= eps;
      for (int j=0;j<nDofs_;j++){
        if (abs(colVal(j)) > 1e-30){
          tripletList.push_back( Eigen::Triplet<scalar_type>( j,i,colVal(j) ) );
        }
      }
    }
    jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  // finite difference apply Jacbian, keep for validation
  void applyJacobian_fd(const state_type &U,
		     const dense_matrix_type & A,
		     scalar_type t,
		     dense_matrix_type &JA) const
  {

    const scalar_type eps = 1.e-5;
    state_type Up(3*nx_*ny_);
    velocity_type V0(3*nx_*ny_);
    velocity_type V_perturb(3*nx_*ny_);

    velocity(U,t,V0);

    for (int k=0; k < A.cols(); k++){
      Up = U + eps*A.col(k);
      velocity(Up,t,V_perturb);
      auto lcol = JA.col(k);
      lcol = 1./eps*(V_perturb - V0 ) ;
    }
  }

protected:
  void setup()
  {
    dx_ = Lx_/static_cast<scalar_type>(nx_);
    dy_ = Ly_/static_cast<scalar_type>(ny_);
    nDofs_ = 3*nx_*ny_;

    const int N_cell = nx_*ny_;
    xGrid_.resize(N_cell);
    yGrid_.resize(N_cell);
    for (int j=0; j < ny_; j++){
      for (int i=0; i < nx_; i ++){
        xGrid_(index_mapper(i,j)) = dx_*i + dx_*0.5;
        yGrid_(index_mapper(i,j)) = dy_*j + dy_*0.5;
      }
    }
    jac_ = createJacobian();
  };

private:
  // mapping from k (state index), i (x-index), and  (y-index) to global index
  int state_index_mapper(const int k,const int i,const int j) const
  {
    const int state_global_indx = 3 * ( (modulus(j,ny_))*nx_ + modulus(i,nx_) ) + k;
    return state_global_indx;
  }

  // mapping from i (x-index), and  (y-index) to global index
  int index_mapper(const int i,const int j) const
  {
    const int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
    return global_indx;
  }

  std::array<int,2> get_ij_from_gid(const int gid) const
  {
      const int j = gid/nx_;
      const int i = gid%nx_;
      std::array<int,2> ij = {i,j};
      return ij;
  }

  int modulus(const int n,const int M) const {
    return ((n % M) + M) % M;
  }


private:
  int nx_ = 64;
  int ny_ = nx_;
  scalar_type Lx_ = static_cast<scalar_type>(5);
  scalar_type Ly_ = Lx_;
  scalar_type dx_ = Lx_/static_cast<scalar_type>(nx_);
  scalar_type dy_ = Ly_/static_cast<scalar_type>(ny_);

  eigVec xGrid_; // mesh points coordinates
  eigVec yGrid_; // mesh points coordinates

  scalar_type g_ = 6.;
  scalar_type mu_ic_ = 0.125;
  scalar_type mu_f_ = 0.2;
  int nDofs_ = {};

  mutable std::vector<Eigen::Triplet<scalar_type>> tripletList = {};
  mutable jacobian_type jac_ = {};
};

}}
#endif
