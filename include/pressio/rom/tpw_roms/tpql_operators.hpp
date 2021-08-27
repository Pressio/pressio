#include "pressio/ops.hpp"

#ifndef TPQL_OPS_HPP_
#define TPQL_OPS_HPP_

namespace pressio{ namespace rom{ namespace experimental{

template <typename scalar_t, typename matrix_t>
Eigen::Matrix<scalar_t, -1 , 1> getCol(const matrix_t & Phi, int i){
  Eigen::Matrix<scalar_t,-1,1> PhiCol = Phi.col(i);
  return PhiCol;
}

template <typename matrix_t, typename vec_t>
void setCol(matrix_t & Phi,const vec_t & PhiCol, int i){
  Phi.col(i) = PhiCol;
}




template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
std::vector< std::vector< std::vector< scalar_t > > > computeBasisTransposeTimesMixedHessian(const app_t & appObj,const state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi)
{
  auto f0 = appObj.createVelocity();
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  scalar_t eps = 1e-3;
  auto yp = ::pressio::ops::clone(y);
  auto muP = ::pressio::ops::clone(mu);

  std::vector< std::vector< std::vector< scalar_t > > > BasisTransposeTimesMixedHessian(romDim);
  Eigen::Matrix<scalar_t,-1,1> BasisTransposeTimesMixedHessianCol(romDim);

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesMixedHessian[i].resize(romDim);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < romDim; j++){
      BasisTransposeTimesMixedHessian[i][j].resize(numParams);
    }
  }
  for  (int i = 0; i < romDim; i++){
    auto PhiCol = getCol<scalar_t>(Phi,i);
    for (int j=0; j < numParams; j++){
      ::pressio::ops::set_zero(workingVector);
      ::pressio::ops::update(yp,0.,y,1.,PhiCol,eps);
      muP(j) += eps;
      //appObj.updateParams(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,0.,ftmp,1.);
      muP(j) -= eps;

      muP(j) -= eps;
      //appObj.updateParams(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      workingVector -= ftmp;
      ::pressio::ops::update(workingVector,1.,ftmp,-1.);
      muP(j) += eps;

      yp = y - eps*PhiCol;
      ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
      muP(j) += eps;
      //appOBj.updateParams(muP)
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.);
      muP(j) -= eps;

      muP(j) -= eps;
      //appObj.udpateParams(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,1.);
      muP(j) += eps;

      ::pressio::ops::scale(workingVector,0.25/(eps*eps));

      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0.,BasisTransposeTimesMixedHessianCol);
      for (int k = 0; k < romDim; k++){
        BasisTransposeTimesMixedHessian[k][i][j] = BasisTransposeTimesMixedHessianCol(k);
      } 
    }
  }
  return BasisTransposeTimesMixedHessian; 
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
std::vector< std::vector< std::vector< scalar_t > > > computeBasisTransposeTimesParameterHessian(const app_t & appObj,const state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi)
{
  auto f0 = appObj.createVelocity();
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  scalar_t eps = 1e-3;

  auto muP = ::pressio::ops::clone(mu);

  std::vector< std::vector< std::vector< scalar_t > > > BasisTransposeTimesParameterHessian(romDim);
  Eigen::Matrix<scalar_t , -1,1> BasisTransposeTimesParameterHessianCol(romDim);
  for (int i=0; i < romDim; i++){
    BasisTransposeTimesParameterHessian[i].resize(numParams);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < numParams; j++){
      BasisTransposeTimesParameterHessian[i][j].resize(numParams);
    }
  }
  for  (int i = 0; i < numParams; i++){
    for (int j=0; j < numParams; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        muP(i) += 2.*eps;
        //appObj.updateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.);
        muP(i) -= 2.*eps;

        muP(i) += eps;
        //appObj.updateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.);
        muP(i) -= eps;

        ::pressio::ops::update(workingVector,1.,f0,-30.);
         
        muP(i) -= eps;
        //appObj.updateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.);
        muP(i) += eps;

        muP(i) -= 2.*eps; 
        //appObj.updateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);
        muP(i) += 2.*eps;
        ::pressio::ops::scale(workingVector, 1./(12.*(eps*eps)));
      }
      else{
        muP(i) += eps;
        muP(j) += eps;
        //appObj.updateParams(muP);
        ftmp.setZero();
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.);
        muP(i) -= eps;
        muP(j) -= eps;

        muP(i) += eps;
        muP(j) -= eps;
        //appObj.updateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);
        muP(i) -= eps;
        muP(j) += eps;

        muP(i) -= eps;
        muP(j) += eps;
        //appOBj.updateParams(muP)
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);
        muP(i) += eps;
        muP(j) -= eps;

        muP(i) -= eps;
        muP(j) -= eps;
        //appObj.udpateParams(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.);
        muP(i) += eps;
        muP(j) += eps;

        ::pressio::ops::scale(workingVector,0.25/(eps*eps));

      }
      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0.,BasisTransposeTimesParameterHessianCol);

      for (int k = 0; k < romDim; k++){
        BasisTransposeTimesParameterHessian[k][i][j] = BasisTransposeTimesParameterHessianCol(k);
      } 
    }
  }
  return BasisTransposeTimesParameterHessian; 
}



template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
Eigen::MatrixXd computeBasisTransposeTimesParameterJacobian(const app_t & appObj,const state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi)
{
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto numParams = ::pressio::ops::extent(mu,0);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  scalar_t eps = 1e-3;
  Eigen::Matrix<scalar_t, -1, -1> PhiTJPhi(romDim,numParams);
  Eigen::Matrix<scalar_t, -1, 1> PhiTJPhiCol(romDim);

  for  (int i = 0; i < numParams; i++){
    ::pressio::ops::set_zero(workingVector);

    Eigen::Matrix<scalar_t, -1, 1>  muP = mu;
    muP(i) += eps;
    //appObj.updateParams(muP); 
    ::pressio::ops::set_zero(ftmp);
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,0.,ftmp,1.);

    muP(i) -= 2.*eps;
    ::pressio::ops::set_zero(ftmp);
    //appObj.updateParams(muP); 
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,1.,ftmp,-1.);

    ::pressio::ops::scale(workingVector, 0.5/eps);

    PhiTJPhi.col(i) = Phi.transpose() * workingVector;
    ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., PhiTJPhiCol);
    setCol(PhiTJPhi,PhiTJPhiCol,i);
  
  }
  return PhiTJPhi;
} 


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
Eigen::Matrix<scalar_t,-1,1> computeBasisTransposeTimesVelocity(const app_t & appObj, const state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f = appObj.createVelocity();
  appObj.velocity(y,t,f);
  Eigen::Matrix<scalar_t, -1, 1> PhiTfCol(romDim);

  ::pressio::ops::product(::pressio::transpose(),1., Phi,f, 0.,PhiTfCol);
  return PhiTfCol;
} 



template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
Eigen::MatrixXd computeBasisTransposeTimesJacobianTimesBasis(const app_t & appObj, const state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto fp = appObj.createVelocity();
  auto yp = ::pressio::ops::clone(y);
  auto workingVector = appObj.createVelocity();

  // Finite difference step size
  scalar_t eps = 1e-3;

  // Create result
  Eigen::Matrix<scalar_t, -1, -1> PhiTJPhi(romDim,romDim);
  Eigen::Matrix<scalar_t, -1, 1> PhiTJPhiCol(romDim);


  // Finite difference loop
  for  (int i = 0; i < romDim; i++){
    ::pressio::ops::set_zero(workingVector);

    auto PhiCol = getCol<scalar_t>(Phi,i);
    ::pressio::ops::update(yp,0.,y,1.,PhiCol,eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,0.,fp,1.);

    ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,1.,fp,-1.);

    ::pressio::ops::scale(workingVector,0.5/eps); 
 
    ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., PhiTJPhiCol);
    setCol(PhiTJPhi,PhiTJPhiCol,i);
  }
  return PhiTJPhi;
} 
  
template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
std::vector< std::vector< std::vector< scalar_t > > > computeBasisTransposeTimesHessianTimesBasisTimesBasis(const app_t & appObj,const state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi)
{
  auto romDim = Phi.cols();
  auto fomDim = y.size();

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  scalar_t eps = 1e-3;
  auto yp = ::pressio::ops::clone(y); //pressio ops clone 

  std::vector< std::vector< std::vector< scalar_t > > > BasisTransposeTimesHessianTimesBasisTimesBasis(romDim);
  Eigen::Matrix<scalar_t, -1, 1> HCol(romDim);

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesHessianTimesBasisTimesBasis[i].resize(romDim);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < romDim; j++){
      BasisTransposeTimesHessianTimesBasisTimesBasis[i][j].resize(romDim);
    }
  }

  for  (int i = 0; i < romDim; i++){
    for (int j=0; j < romDim; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        auto PhiCol = getCol<scalar_t>(Phi,i);

        //yp = y + 2.*eps*PhiCol;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol,2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.);

        //yp = y + eps*PhiCol (or yp = yp - eps PhiCol);
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.);

    
        ::pressio::ops::update(workingVector,1.,f0,-30.);
    
      
        //yp = y - eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.);
       
        //yp = y - 2.*eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);
 
        ::pressio::ops::scale(workingVector,1./(12.*(eps*eps)));
      }
      else{
        auto PhiCol_i = getCol<scalar_t>(Phi,i);
        auto PhiCol_j = getCol<scalar_t>(Phi,j);

        //yp = y + eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps,PhiCol_j,eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.);

        //yp = y + eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);

        //yp = y - eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_i,-2.*eps,PhiCol_j,2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.);
       
        //yp = y - eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.);
 
        ::pressio::ops::scale(workingVector,0.25/(eps*eps));

      }
      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., HCol);
      for (int k = 0; k < romDim; k++){
        BasisTransposeTimesHessianTimesBasisTimesBasis[k][i][j] = HCol(k);
      } 
    }
  }
  return BasisTransposeTimesHessianTimesBasisTimesBasis; 
} 

}}}
#endif
