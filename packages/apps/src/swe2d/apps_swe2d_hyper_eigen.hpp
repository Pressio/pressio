
#ifndef PRESSIO_APPS_SWE_HYPER_HPP_
#define PRESSIO_APPS_SWE_HYPER_HPP_

#include "Eigen/Core"
#include "rusanovflux.hpp"
#include <unordered_map>
#include <chrono>

namespace pressio{ namespace apps{

template <typename gids_t>
class swe2d_hyper
{
  using eigVec = Eigen::VectorXd;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using eig_sp_mat = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using jacobian_type	= eig_sp_mat;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  swe2d_hyper(const scalar_type Lx,
	      const scalar_type Ly,
              const int nx,
	      const int ny,
	      const scalar_type params[3],
              const gids_t & gidsSample,
	      const gids_t & gidsStencil)
    : nx_(nx), ny_(ny),
      Lx_(Lx), Ly_(Ly),
      g_(params[0]),
      mu_ic_(params[1]),
      mu_f_(params[2]),
      gidsSample_(gidsSample),
      gidsStencil_(gidsStencil)
  {
    this->setup();
  }

  swe2d_hyper(const int n,
	      const scalar_type params[3],
              const gids_t gidsSample,
	      const gids_t gidsStencil)
    : nx_(n), ny_(n),
      g_(params[0]),
      mu_ic_(params[1]),
      mu_f_(params[2]),
      gidsSample_(gidsSample),
      gidsStencil_(gidsStencil)
  {
    this->setup();
  }

  // Get Gaussian IC at full mesh
  state_type getGaussianICFull(const scalar_type mu) const
  {
    state_type U0(3*nx_*ny_);
    for (int j=0; j < ny_; j++){
      for (int i=0; i < nx_; i++){
        const auto sid = state_index_mapper(0,i,j);
        U0(sid) = 1. + mu_ic_*exp( - ( pow( xGrid_(index_mapper(i,j)) - 1.5,2)
				       + pow( yGrid_(index_mapper(i,j)) - 1.5,2) ));
        U0(sid+1) = 0.;
        U0(sid+2) = 0.;
      }
    }
    return U0;
  }

  // Get Gaussian IC on stencil mesh
  state_type getGaussianIC(const scalar_type mu) const
  {
    state_type U0(3*nDofsStencil_);
    for (int smpsIndexCounter=0; smpsIndexCounter < nDofsStencil_; smpsIndexCounter++)
      {
        const auto smpsGlobalIndex = smpsLidToGidMap_(smpsIndexCounter);
        const auto ij = get_ij_from_gid(smpsGlobalIndex);
        const auto i = ij[0];
        const auto j = ij[1];
        const auto indexStencil = 3*smpsIndexCounter;
        U0(indexStencil) = 1. + mu_ic_*exp( - ( pow( xGrid_(index_mapper(i,j)) - 1.5,2)
                                   + pow( yGrid_(index_mapper(i,j)) - 1.5,2) ));
        U0(indexStencil+1) = 0.;
        U0(indexStencil+2) = 0.;
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
    velocity_type V(3*nDofsSample_);
    return V;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3*nDofsSample_,3*nDofsStencil_);
    return J;
  }

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & A) const
  {
    jacobian_type JA(3*nDofsSample_, A.cols() );
    return JA;
  }

  void velocity(const state_type & U,
                const scalar_type & /*t*/,
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
    for (int smIndexCounter=0; smIndexCounter < nDofsSample_; smIndexCounter++)
      {
        const auto sampleMeshIndex = smLidToGidMap_(smIndexCounter);
        const auto ij = get_ij_from_gid(sampleMeshIndex); // get ij location
        const auto i = ij[0];
        const auto j = ij[1];

        isL = 3*(smpsGidToLidMap_.find(index_mapper(i-1,j))->second);
	isR = 3*(smpsGidToLidMap_.find(index_mapper(i,j))->second);
	rusanovFluxFullStateIn(FL, U, isL, isR, nx, g_);

	isL = 3*(smpsGidToLidMap_.find(index_mapper(i,j))->second);
	isR = 3*(smpsGidToLidMap_.find(index_mapper(i+1,j))->second);
	rusanovFluxFullStateIn(FR, U, isL, isR, nx, g_);

	isD = 3*(smpsGidToLidMap_.find(index_mapper(i,j-1))->second);
	isU = 3*(smpsGidToLidMap_.find(index_mapper(i,j))->second);
	rusanovFluxFullStateIn(FD, U, isD, isU, ny, g_);

	isD = 3*(smpsGidToLidMap_.find(index_mapper(i,j))->second);
	isU = 3*(smpsGidToLidMap_.find(index_mapper(i,j+1))->second);
	rusanovFluxFullStateIn(FU, U, isD, isU, ny, g_);

        const auto isSample  = 3*(smGidToLidMap_.find(index_mapper(i,j))->second);
        const auto isStencil = 3*(smpsGidToLidMap_.find(index_mapper(i,j))->second);

        forcing[0] = 0;
        forcing[1] = mu_f_*U(isStencil+2);
        forcing[2] = mu_f_*U(isStencil+1);
	for (int k=0;k<3;k++){
	  V(isSample+k) = -1./dx_*(FR[k] - FL[k]) - 1./dy_*(FU[k] - FD[k]) + forcing[k];
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

      const std::array<scalar_type, 2> nx = { 1, 0 };
      const std::array<scalar_type, 2> ny = { 0, 1 };
      tripletList.clear();

      for (int smIndexCounter=0; smIndexCounter < nDofsSample_; smIndexCounter++)
      {
	const auto sampleMeshGlobalIndex = smLidToGidMap_(smIndexCounter);
	const auto ij = get_ij_from_gid(sampleMeshGlobalIndex);
	const auto smps_lid =     3*smpsGidToLidMap_.find(index_mapper(ij[0]    ,ij[1]    ))->second;
	const auto smps_lid_im1 = 3*smpsGidToLidMap_.find(index_mapper(ij[0] - 1,ij[1]    ))->second;
	const auto smps_lid_ip1 = 3*smpsGidToLidMap_.find(index_mapper(ij[0] + 1,ij[1]    ))->second;
	const auto smps_lid_jm1 = 3*smpsGidToLidMap_.find(index_mapper(ij[0]    ,ij[1] - 1))->second;
	const auto smps_lid_jp1 = 3*smpsGidToLidMap_.find(index_mapper(ij[0]    ,ij[1] + 1))->second;

	rusanovFluxJacobianFullStateIn(JL_L,JR_L,U,smps_lid_im1, smps_lid    ,nx,g_);
	rusanovFluxJacobianFullStateIn(JL_R,JR_R,U,smps_lid    , smps_lid_ip1,nx,g_);
	rusanovFluxJacobianFullStateIn(JD_D,JU_D,U,smps_lid_jm1, smps_lid    ,ny,g_);
	rusanovFluxJacobianFullStateIn(JD_U,JU_U,U,smps_lid    , smps_lid_jp1,ny,g_);

	const auto sm_lid = smIndexCounter*3;
	for (int j = 0; j < 3; ++j) {
	  for (int i = 0; i < 3; ++i) {
	    auto val = -1./dx_*(JL_R[i][j] - JR_L[i][j]) - 1./dy_*(JD_U[i][j] - JU_D[i][j]);
	    tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+i,smps_lid+j,val ) );

	    val =  1./dx_*JL_L[i][j];
	    tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+i,smps_lid_im1+j,val ) );

	    val = -1./dx_*JR_R[i][j];
	    tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+i,smps_lid_ip1+j,val ) );

	    val = 1./dy_*JD_D[i][j];
	    tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+i,smps_lid_jm1+j,val ) );

	    val = -1./dy_*JU_U[i][j];
	    tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+i,smps_lid_jp1+j,val ) );

	  }
	}

        // add forcing terms
        auto val = mu_f_;
        tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+1,smps_lid+2,val ) );
        tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid+2,smps_lid+1,val ) );
      }
      jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  // Finite difference version of Jacobian, keep for validation
  void jacobian_fd(const state_type & U,
                  const scalar_type /*t*/,
                  jacobian_type & jac) const
  {
    velocity_type V0(3*nDofsSample_);
    velocity_type Vp(3*nDofsSample_);
    state_type Up(U);
    velocity(U,0.,V0);
    const scalar_type eps = 1e-5;
    tripletList.clear();
    // Loop through all points on stencil mesh and perturb
    for (int i = 0; i < nDofsStencil_; i++){
      for (int k=0; k<3; k++){
        auto smps_lid = 3*i + k;
        Up(smps_lid) += eps;
        velocity(Up,0.,Vp);
        auto colVal = (Vp - V0)/eps;
        Up(smps_lid) -= eps;
        // Loop through all dofs on sample mesh and add to triplet if non-zero
        for (int j=0;j<nDofsSample_;j++){
          for (int l=0;l<3;l++){
            // Get index to use for velocity
            auto sm_lid = 3*j + l;
            if (abs(colVal(sm_lid)) > 1e-30){
              tripletList.push_back( Eigen::Triplet<scalar_type>( sm_lid,smps_lid,colVal(sm_lid) ) );
            }
          }
        }
      }
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

  void applyJacobian_fd(const state_type &U,
		     const dense_matrix_type & A,
		     const scalar_type t,
		     dense_matrix_type &JA) const
  {
    const scalar_type eps = 1.e-5;
    state_type Up(3*nDofsStencil_);
    velocity_type V0(3*nDofsSample_);
    velocity_type V_perturb(3*nDofsSample_);

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
    nDofsSample_ = gidsSample_.size();
    smLidToGidMap_.resize(nDofsSample_);
    for (int i = 0; i < nDofsSample_; i++)
    {
      smLidToGidMap_(i) = gidsSample_[i];
      smGidToLidMap_[gidsSample_[i]] = i;
    }

    nDofsStencil_ = gidsStencil_.size();
    smpsLidToGidMap_.resize(nDofsStencil_);
    for (int i = 0; i < nDofsStencil_; i++)
    {
      smpsLidToGidMap_(i) = gidsStencil_[i];
      smpsGidToLidMap_[gidsStencil_[i]] = i;
    }

    dx_ = Lx_/static_cast<scalar_type>(nx_);
    dy_ = Ly_/static_cast<scalar_type>(ny_);
    int N_cell = nx_*ny_;
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

  // mapping from k (state index), i (x-index), and  (y-index) to global index
  int state_index_mapper(const int k,const int i,const int j) const {
    const int state_global_indx = 3 * ( (modulus(j,ny_))*nx_ + modulus(i,nx_) ) + k;
    return state_global_indx;
  }

  // mapping from i (x-index), and  (y-index) to global index
  int index_mapper(const int i,const int j) const {
    const int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
    return global_indx;
  }

  // get ij index from global index
  std::array<int,2> get_ij_from_gid(const int gid) const {
      const int j = gid/nx_;
      const int i = gid%nx_;
      const std::array<int,2> ij = {i,j};
      return ij;
  }

  int modulus(int n,int M) const {
    return ((n % M) + M) % M;
  }

private:
  int nx_ = {};
  int ny_ = {};
  scalar_type Lx_ = 5.;
  scalar_type Ly_ = Lx_;
  scalar_type dx_ = {};
  scalar_type dy_ = {};
  int nDofsSample_ = {};
  int nDofsStencil_ = {};

  eigVec xGrid_;
  eigVec yGrid_;

  scalar_type g_ = {};
  scalar_type mu_ic_ = {};
  scalar_type mu_f_ = {};
  gids_t gidsSample_ = {};
  gids_t gidsStencil_ = {};

  Eigen::Matrix<int,-1,1> smLidToGidMap_;
  Eigen::Matrix<int,-1,1> smpsLidToGidMap_;

  using mapping_t = std::unordered_map<int,int>;
  mapping_t smGidToLidMap_;
  mapping_t smpsGidToLidMap_;

  mutable std::vector<  Eigen::Triplet<scalar_type> > tripletList;
  mutable jacobian_type jac_;
};

}}
#endif
