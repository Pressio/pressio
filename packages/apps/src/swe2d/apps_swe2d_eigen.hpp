
#ifndef SWE_HPP_
#define SWE_HPP_
#include "Eigen/Core"
#include "rusanovflux.hpp"
namespace pressio{ namespace apps{
class swe2d 
{

  using eigVec = Eigen::VectorXd;
  using ui_t = unsigned int;
protected:

  eigVec U_;
  eigVec xGrid_; // mesh points coordinates
  eigVec yGrid_; // mesh points coordinates

public:
  using scalar_type	= double; 
  using scalar_t        = scalar_type;
  using state_type	= eigVec;
  using state_t          = state_type;
  using velocity_type	= eigVec;
  using eig_sp_mat = Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, int>;
  using jacobian_type	= eig_sp_mat;
  using jacobian_t = jacobian_type;
  using velocity_t = velocity_type;
  using dense_matrix_type = Eigen::MatrixXd; 


private:
  mutable std::vector<  Eigen::Triplet<scalar_t> > tripletList;
  scalar_t Lx_;
  scalar_t Ly_;
  int nx_;
  int ny_;
  scalar_t dx_;
  scalar_t dy_;
  int nDofs_;
  scalar_t g_;
  scalar_t mu_ic_;
  scalar_t mu_f_;
  mutable jacobian_t jac_;

  // mapping from k (state index), i (x-index), and  (y-index) to global index
  int state_index_mapper(int k, int i,int j) const {
    int state_global_indx = 3 * ( (modulus(j,ny_))*nx_ + modulus(i,nx_) ) + k;
    return state_global_indx;
  }

  // mapping from i (x-index), and  (y-index) to global index
  int index_mapper(int i,int j) const {
    int global_indx = (modulus(j,ny_))*nx_ + modulus(i,nx_);
    return global_indx;
  }

  std::array<int,2> get_ij_from_gid(int gid) const {
      int j = gid/nx_;
      int i = gid%nx_;
      std::array<int,2> ij = {i,j};
      return ij;
  }

  int modulus(int n,int M) const {
    return ((n % M) + M) % M;
  }


public:
  swe2d(scalar_t Lx, scalar_t Ly, int nx, int ny, scalar_t params[3]) :
    nx_(nx),Lx_(Lx),ny_(ny),Ly_(Ly), g_(params[0]),mu_ic_(params[1]),mu_f_(params[2]){
      this->setup();
    }


  //==========
  void velocity(const state_type & U, const scalar_t & t, velocity_type & V) const
  {
    scalar_t FL[3];
    scalar_t FR[3];
    scalar_t FU[3];
    scalar_t FD[3];
    scalar_t forcing[3];
    std::array<scalar_t, 2> nx = { 1, 0 };
    std::array<scalar_t, 2> ny = { 0, 1 };

    int isL;
    int isR;
    int isU;
    int isD;
    for (int i=0; i < nx_  ; i++){
      for (int j=0; j < ny_ ; j++){
        isL = state_index_mapper(0,i-1,j);
        isR = state_index_mapper(0,i,j);
	rusanovflux_eigen(FL,U,isL,isR,nx,g_);
        isL = state_index_mapper(0,i,j);
        isR = state_index_mapper(0,i+1,j);
	rusanovflux_eigen(FR,U,isL,isR,nx,g_);

        isD = state_index_mapper(0,i,j-1);
        isU = state_index_mapper(0,i,j);
	rusanovflux_eigen(FD,U,isD,isU,ny,g_);

        isD = state_index_mapper(0,i,j);
        isU = state_index_mapper(0,i,j+1);
	rusanovflux_eigen(FU,U,isD,isU,ny,g_);

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
  /*

   */
  void jacobian(const state_type & U, const scalar_type /*t*/, jacobian_type & jac) const{
      scalar_t *UL;
      scalar_t *UR;
      scalar_t *UU;
      scalar_t *UD;
      scalar_t *V_view;

      scalar_t JL_L[3][3];
      scalar_t JR_L[3][3];
      scalar_t JL_R[3][3];
      scalar_t JR_R[3][3];
      scalar_t JD_D[3][3];
      scalar_t JU_D[3][3];
      scalar_t JD_U[3][3];
      scalar_t JU_U[3][3];

      std::array<scalar_t, 2> nx = { 1, 0 };
      std::array<scalar_t, 2> ny = { 0, 1 };


      //=================
      tripletList.clear();
      for (int sid=0; sid < nx_*ny_  ; sid++){
        auto ij = get_ij_from_gid(sid);
        auto gid = state_index_mapper(0,ij[0],ij[1]);
        auto gid_im1 = state_index_mapper(0,ij[0] - 1,ij[1]);
        auto gid_ip1 = state_index_mapper(0,ij[0] + 1,ij[1]);
        auto gid_jm1 = state_index_mapper(0,ij[0],ij[1] - 1);
        auto gid_jp1 = state_index_mapper(0,ij[0],ij[1] + 1);

        rusanovflux_jacobian_eigen(JL_L,JR_L,U,gid_im1,gid,nx,g_);
        rusanovflux_jacobian_eigen(JL_R,JR_R,U,gid,gid_ip1,nx,g_);
        rusanovflux_jacobian_eigen(JD_D,JU_D,U,gid_jm1,gid,ny,g_);
        rusanovflux_jacobian_eigen(JD_U,JU_U,U,gid,gid_jp1,ny,g_);
        for (int j = 0; j < 3; ++j) {
          for (int i = 0; i < 3; ++i) {
            auto val = -1./dx_*(JL_R[i][j] - JR_L[i][j]) - 1./dy_*(JD_U[i][j] - JU_D[i][j]);
            tripletList.push_back( Eigen::Triplet<scalar_t>( gid+i,gid+j,val ) );

            val =  1./dx_*JL_L[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_t>( gid+i,gid_im1+j,val ) );

            val = -1./dx_*JR_R[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_t>( gid+i,gid_ip1+j,val ) );

            val = 1./dy_*JD_D[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_t>( gid+i,gid_jm1+j,val ) );

            val = -1./dy_*JU_U[i][j];
            tripletList.push_back( Eigen::Triplet<scalar_t>( gid+i,gid_jp1+j,val ) );

          }
        }
        // add forcing terms
        auto val = mu_f_;
        tripletList.push_back( Eigen::Triplet<scalar_t>( gid+1,gid+2,val ) );
        tripletList.push_back( Eigen::Triplet<scalar_t>( gid+2,gid+1,val ) );
      }
      jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void jacobian_fd(const state_type & U, const scalar_type /*t*/, jacobian_type & jac) const{
    velocity_type V0(nDofs_);
    velocity_type Vp(nDofs_);
    state_type Up(U);
    velocity(U,0.,V0);
    scalar_t eps = 1e-5;
    tripletList.clear();
    for (int i = 0; i < nDofs_; i++){
      Up(i) += eps;
      velocity(Up,0.,Vp);
      auto colVal = (Vp - V0)/eps;
      Up(i) -= eps;
      for (int j=0;j<nDofs_;j++){
        if (abs(colVal(j)) > 1e-30){
          tripletList.push_back( Eigen::Triplet<scalar_t>( j,i,colVal(j) ) );
        }
      } 
    }
    jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }


  void applyJacobian(const state_type &U,
		     const dense_matrix_type & A,
		     scalar_t t,
		     dense_matrix_type &JA) const
  {
    jacobian(U,t,jac_);
    JA = jac_*A;
  }

  void applyJacobian_fd(const state_type &U,
		     const dense_matrix_type & A,
		     scalar_t t,
		     dense_matrix_type &JA) const
  {
    
    scalar_t eps = 1.e-5;
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

  //========
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

protected:
  void setup(){
    // distribute cells
    dx_ = Lx_/static_cast<scalar_t>(nx_);
    dy_ = Ly_/static_cast<scalar_t>(ny_);
    nDofs_ = 3*nx_*ny_;
    int N_cell = nx_*ny_;
    xGrid_.resize(N_cell);
    yGrid_.resize(N_cell);
    for (int i=0; i < nx_; i ++){
      for (int j=0; j < ny_; j++){
        xGrid_(index_mapper(i,j)) = dx_*i + dx_*0.5;
        yGrid_(index_mapper(i,j)) = dy_*j + dy_*0.5;
      }
    }
    U_.resize(N_cell*3);
    jac_ = createJacobian();
  };

public:
  state_type getGaussianIC(scalar_t mu) const{
    state_type U0(3*nx_*ny_);
    for (int i=0; i < nx_; i++){
      for (int j=0; j < ny_; j++){
        auto sid = state_index_mapper(0,i,j);
        U0(sid) = 1. + mu_ic_*exp( - ( pow( xGrid_(index_mapper(i,j)) - 1.5,2)
                                   + pow( yGrid_(index_mapper(i,j)) - 1.5,2) ));
        U0(sid+1) = 0.;
        U0(sid+2) = 0.;
      }
    }
    return U0;
  }
};
}}
#endif
