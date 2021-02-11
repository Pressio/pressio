
#ifndef SWE_HYPER_HPP_
#define SWE_HYPER_HPP_
#include "Eigen/Core"
#include "rusanovflux.hpp"
namespace pressio{ namespace apps{
template <typename gids_t>
class swe2d_hyper
{

  using eigVec = Eigen::VectorXd;
protected:

  eigVec xGrid_; // mesh points coordinates
  eigVec yGrid_; // mesh points coordinates

public:
  using scalar_type	= double; 
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using eig_sp_mat = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;
  using jacobian_type	= eig_sp_mat;
  using dense_matrix_type = Eigen::MatrixXd; 

private:
  mutable std::vector<  Eigen::Triplet<scalar_type> > tripletList;
  scalar_type Lx_;
  scalar_type Ly_;
  int nx_;
  int ny_;
  scalar_type dx_;
  scalar_type dy_;
  int nDofsSample_;
  int nDofsStencil_;

  mutable scalar_type g_;
  mutable scalar_type mu_ic_;
  mutable scalar_type mu_f_;
  gids_t gidsSample_;
  gids_t gidsStencil_;
  mutable jacobian_type jac_;
  mutable std::map<int,int> smLidToGidMap_;
  mutable std::map<int,int> smGidToLidMap_;
  mutable std::map<int,int> smpsLidToGidMap_;
  mutable std::map<int,int> smpsGidToLidMap_;


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


public:
  swe2d_hyper(const scalar_type Lx,const scalar_type Ly, 
              const int nx,const int ny,const scalar_type params[3], 
              const gids_t gidsSample,const gids_t gidsStencil)
    : nx_(nx),Lx_(Lx),ny_(ny),Ly_(Ly), 
      g_(params[0]),mu_ic_(params[1]),
      mu_f_(params[2]), gidsSample_(gidsSample),
      gidsStencil_(gidsStencil)
    {
        this->setup();
    }

  // Sets parameter values
  void setParams(const scalar_type params[3]) const
  {
    g_ = params[0];
    mu_ic_ = params[1]; 
    mu_f_ = params[2];
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
    for (int sampleMeshIndexCounter=0; sampleMeshIndexCounter < nDofsSample_  ; sampleMeshIndexCounter++){
        const auto sampleMeshIndex = smLidToGidMap_[sampleMeshIndexCounter];
        const auto ij = get_ij_from_gid(sampleMeshIndex); // get ij location
        const auto i = ij[0];
        const auto j = ij[1];
        isL = 3*smpsGidToLidMap_[index_mapper(i-1,j)]; //gives global id
        isR = 3*smpsGidToLidMap_[index_mapper(i,j)];

	rusanovFluxFullStateIn(FL,U,isL,isR,nx,g_);
        isL = 3*smpsGidToLidMap_[index_mapper(i,j)];
        isR = 3*smpsGidToLidMap_[index_mapper(i+1,j)];
	rusanovFluxFullStateIn(FR,U,isL,isR,nx,g_);

        isD = 3*smpsGidToLidMap_[index_mapper(i,j-1)];
        isU = 3*smpsGidToLidMap_[index_mapper(i,j)];
	rusanovFluxFullStateIn(FD,U,isD,isU,ny,g_);

        isD = 3*smpsGidToLidMap_[index_mapper(i,j)];
        isU = 3*smpsGidToLidMap_[index_mapper(i,j+1)];
	rusanovFluxFullStateIn(FU,U,isD,isU,ny,g_);

        const auto isSample = 3*smGidToLidMap_[index_mapper(i,j)];
        const auto isStencil = 3*smpsGidToLidMap_[index_mapper(i,j)];

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
      for (int sampleMeshIndexCounter=0; sampleMeshIndexCounter < nDofsSample_  ; sampleMeshIndexCounter++){
        const auto sampleMeshGlobalIndex = smLidToGidMap_[sampleMeshIndexCounter];  
        const auto ij = get_ij_from_gid(sampleMeshGlobalIndex);
        const auto smps_lid =     3*smpsGidToLidMap_[index_mapper(ij[0],ij[1])];
        const auto smps_lid_im1 = 3*smpsGidToLidMap_[index_mapper(ij[0] - 1,ij[1])];
        const auto smps_lid_ip1 = 3*smpsGidToLidMap_[index_mapper(ij[0] + 1,ij[1])];
        const auto smps_lid_jm1 = 3*smpsGidToLidMap_[index_mapper(ij[0],ij[1] - 1)];
        const auto smps_lid_jp1 = 3*smpsGidToLidMap_[index_mapper(ij[0],ij[1] + 1)];

        rusanovFluxJacobianFullStateIn(JL_L,JR_L,U,smps_lid_im1, smps_lid    ,nx,g_);
        rusanovFluxJacobianFullStateIn(JL_R,JR_R,U,smps_lid    , smps_lid_ip1,nx,g_);
        rusanovFluxJacobianFullStateIn(JD_D,JU_D,U,smps_lid_jm1, smps_lid    ,ny,g_);
        rusanovFluxJacobianFullStateIn(JD_U,JU_U,U,smps_lid    , smps_lid_jp1,ny,g_);
        const auto sm_lid = sampleMeshIndexCounter*3;
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

  // Finite difference version, keep for validation
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

  //========
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

protected:
  void setup(){
    // Create maps for hyperreduction 
    nDofsSample_ = gidsSample_.size();
    nDofsStencil_ = gidsStencil_.size();
    for (int i = 0; i < nDofsSample_; i++){
      smLidToGidMap_.emplace(i,gidsSample_[i]);
      smGidToLidMap_.emplace(gidsSample_[i],i);
    }
    for (int i = 0; i < nDofsStencil_; i++){
      smpsLidToGidMap_.emplace(i,gidsStencil_[i]);
      smpsGidToLidMap_.emplace(gidsStencil_[i],i);
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

public:
  // Get Gaussian IC at full mesh
  state_type getGaussianICFull(const scalar_type mu) const{
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
  state_type getGaussianIC(const scalar_type mu) const{
    state_type U0(3*nDofsStencil_);
    for (int smpsIndexCounter=0; smpsIndexCounter < nDofsStencil_  ; smpsIndexCounter++){
        const auto smpsGlobalIndex = smpsLidToGidMap_[smpsIndexCounter];
        const auto ij = get_ij_from_gid(smpsGlobalIndex); // get ij location
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
};
}}
#endif
