#include "pressio/ops.hpp"
#include "tpql_rom_datatypes.hpp"
#ifndef TPQL_OPS_HPP_
#define TPQL_OPS_HPP_

namespace pressio{ namespace rom{ namespace experimental{



template <typename scalar_t>
auto getCol(const Tpetra::BlockMultiVector<> & Phi,const int i){

  using map_t   = Tpetra::Map<>;

  using go_t    = typename map_t::global_ordinal_type;
  using lo_t    = typename map_t::local_ordinal_type;

  auto meshMap = Phi.getPointMap();
  auto PhiMultiVector = Phi.getMultiVectorView(); 
  auto PhiCol = PhiMultiVector.getVector(i);
  Tpetra::BlockVector<> PhiColBlockVector(*PhiCol , meshMap, 1);
  return PhiColBlockVector;
}


template <typename scalar_t>
auto getCol(const Tpetra::MultiVector<> & Phi,const int i){
  auto PhiCol = Phi.getVector(i);
  return *PhiCol;
}



template <typename scalar_t, typename matrix_t>
Eigen::Matrix<scalar_t, -1 , 1> getCol(const matrix_t & Phi,const int i){
  Eigen::Matrix<scalar_t,-1,1> PhiCol = Phi.col(i);
  return PhiCol;
}

template <typename scalar_t>
void setCol(Eigen::Matrix<scalar_t,-1,-1> & Phi,const Eigen::Matrix<scalar_t,-1,1>  & PhiCol,const int i){
  Phi.col(i) = PhiCol;
}


template <typename scalar_t, typename layout_t, typename exec_space_t>
void setCol(Kokkos::View<scalar_t**,layout_t,exec_space_t> & Phi,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & PhiCol, const int i){
 auto numRows = ::pressio::ops::extent(Phi,0);
 for (int j = 0; j < numRows; j++){
   Phi(j,i) = PhiCol(j);
 }
}


template <typename scalar_t>
void setHessianCol(std::vector<std::vector<std::vector<scalar_t> > > & H , const Eigen::Matrix<scalar_t,-1,1>  & HCol,const int i, const int j )
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H[k][i][j] = HCol(k);
  } 
}

template <typename scalar_t, typename layout_t, typename exec_space_t>
void setHessianCol(Kokkos::View<scalar_t***,layout_t,exec_space_t> & H,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & HCol,const int i,const int j)
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H(k,i,j) = HCol(k);
  } 
}

template <typename scalar_t>
void setSymmetricHessianCol(std::vector<std::vector<std::vector<scalar_t> > > & H , const Eigen::Matrix<scalar_t,-1,1>  & HCol,const int i, const int j )
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H[k][i][j] = HCol(k);
    H[k][j][i] = HCol(k);
  } 
}

template <typename scalar_t, typename layout_t, typename exec_space_t>
void setSymmetricHessianCol(Kokkos::View<scalar_t***,layout_t,exec_space_t> & H,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & HCol,const int i,const int j)
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H(k,i,j) = HCol(k);
    H(k,j,i) = HCol(k);
  } 
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename dense_matrix_t, typename vector_t>
void computeBasisTransposeTimesMixedHessian_calc(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, dense_matrix_t & BasisTransposeTimesMixedHessian, vector_t & BasisTransposeTimesMixedHessianCol, const scalar_t epsilon1, const scalar_t epsilon2)

{
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto y0_norm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  auto yp = ::pressio::ops::clone(y);
  auto muP = ::pressio::ops::clone(mu);

  //scalar_t epsState = 1e-3;
  scalar_t epsParams = epsilon2;

  //scalar_t scale = 0.25/(epsState*epsParams);

  for  (int i = 0; i < romDim; i++){
    auto PhiCol = getCol<scalar_t>(Phi,i);
    auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
    scalar_t epsState = std::sqrt( (1. + y0_norm)*epsilon1)/PhiColNorm;
    scalar_t scale = 0.25 / (epsState * epsParams);
    for (int j=0; j < numParams; j++){
      ::pressio::ops::set_zero(workingVector);
      ::pressio::ops::update(yp,0.,y,1.,PhiCol,epsState);
      muP(j) += epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);
      muP(j) -= epsParams;

      muP(j) -= epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
      muP(j) += epsParams;

      ::pressio::ops::update(yp,1.,PhiCol,-2.*epsState);
      muP(j) += epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
      muP(j) -= epsParams;

      muP(j) -= epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);
      muP(j) += epsParams;

      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0.,BasisTransposeTimesMixedHessianCol);
      setHessianCol(BasisTransposeTimesMixedHessian,BasisTransposeTimesMixedHessianCol,i,j);
    }
  }
  appObj.updateScalarParameters(mu);
}




template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesMixedHessianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi,const scalar_t epsilon1, const scalar_t epsilon2)

{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesMixedHessian("PhiTH_params", romDim,romDim,numParams);
  vector_t HCol("HCol",romDim);
  computeBasisTransposeTimesMixedHessian_calc(appObj,y,mu,t,Phi,BasisTransposeTimesMixedHessian,HCol,epsilon1, epsilon2);
  return BasisTransposeTimesMixedHessian; 
}





template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
std::vector< std::vector< std::vector< scalar_t > > > computeBasisTransposeTimesMixedHessianEigen(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const scalar_t epsilon1, const scalar_t epsilon2)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t = typename rom_data_t::vector_t;
  using dense_hessian_t = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesMixedHessian(romDim);
  vector_t BasisTransposeTimesMixedHessianCol(romDim);
  BasisTransposeTimesMixedHessianCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesMixedHessian[i].resize(romDim);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < romDim; j++){
      BasisTransposeTimesMixedHessian[i][j].resize(numParams);
      for (int k=0; k < numParams; k++){
        BasisTransposeTimesMixedHessian[i][j][k] = 0.;
      }
    }
  }

  computeBasisTransposeTimesMixedHessian_calc(appObj,y,mu,t,Phi,BasisTransposeTimesMixedHessian,BasisTransposeTimesMixedHessianCol,epsilon1,epsilon2);
  return BasisTransposeTimesMixedHessian; 
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeBasisTransposeTimesParameterHessian_calc(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, hessian_t & BasisTransposeTimesParameterHessian, vector_t & BasisTransposeTimesParameterHessianCol, const scalar_t epsilon)
{
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  scalar_t eps = epsilon;

  auto muP = ::pressio::ops::clone(mu);

  scalar_t scale1 = 1./(12.*eps*eps);
  scalar_t scale2 = 0.25/(eps*eps);

  for  (int i = 0; i < numParams; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        muP(i) += 2.*eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.*scale1);
        muP(i) -= 2.*eps;

        muP(i) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale1);
        muP(i) -= eps;

        ::pressio::ops::update(workingVector,1.,f0,-30.*scale1);
         
        muP(i) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale1);
        muP(i) += eps;

        muP(i) -= 2.*eps; 
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale1);
        muP(i) += 2.*eps;
      }
      else{
        muP(i) += eps;
        muP(j) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.*scale2);
        muP(i) -= eps;
        muP(j) -= eps;

        muP(i) += eps;
        muP(j) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale2);
        muP(i) -= eps;
        muP(j) += eps;

        muP(i) -= eps;
        muP(j) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale2);
        muP(i) += eps;
        muP(j) -= eps;

        muP(i) -= eps;
        muP(j) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.*scale2);
        muP(i) += eps;
        muP(j) += eps;
      }
      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0.,BasisTransposeTimesParameterHessianCol);


      setSymmetricHessianCol(BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,i,j);
    }
  }
  //reset parameters
  appObj.updateScalarParameters(mu);
}



template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesParameterHessianEigen(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const scalar_t epsilon)
{

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);
  using vector_t = typename rom_data_t::vector_t;
  using dense_hessian_t = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesParameterHessian(romDim);
  vector_t BasisTransposeTimesParameterHessianCol(romDim);
  BasisTransposeTimesParameterHessianCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesParameterHessian[i].resize(numParams);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < numParams; j++){
      BasisTransposeTimesParameterHessian[i][j].resize(numParams);
      for (int k=0; k < numParams; k++){
        BasisTransposeTimesParameterHessian[i][j][k] = 0.;
      }
    }
  }
  computeBasisTransposeTimesParameterHessian_calc(appObj,y,mu,t,Phi,BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,epsilon);
  return BasisTransposeTimesParameterHessian; 
}

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesParameterHessianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const scalar_t epsilon)
{

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesParameterHessian("PhiTH_params", romDim,numParams,numParams);
  vector_t BasisTransposeTimesParameterHessianCol("HCol",romDim);
  computeBasisTransposeTimesParameterHessian_calc(appObj,y,mu,t,Phi,BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,epsilon);
  return BasisTransposeTimesParameterHessian; 
}


/*--------------------------------------------------
Phi^T J_params computation
--------------------------------------------------*/
template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename dense_matrix_t, typename vector_t>
void computeBasisTransposeTimesParameterJacobian_calc(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, dense_matrix_t & PhiTJPhi, vector_t & PhiTJPhiCol, const scalar_t epsilon)
{
  auto ftmp = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  auto workingVector = appObj.createVelocity();

  auto numParams = ::pressio::ops::extent(mu,0);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  scalar_t eps = epsilon;
  scalar_t scale = 0.5/eps;
  auto muP = ::pressio::ops::clone(mu);
  for  (int i = 0; i < numParams; i++){
    ::pressio::ops::set_zero(workingVector);

    muP(i) += eps;
    appObj.updateScalarParameters(muP); 
    ::pressio::ops::set_zero(ftmp);
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);

    muP(i) -= 2.*eps;
    ::pressio::ops::set_zero(ftmp);
    appObj.updateScalarParameters(muP); 
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

    ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., PhiTJPhiCol);
    setCol(PhiTJPhi,PhiTJPhiCol,i);
    std::cout << "PhiTJPhiCol norm = " << ::pressio::ops::norm2(PhiTJPhiCol) << std::endl;
  
  }
  // reset parameters to original state
  appObj.updateScalarParameters(mu);
} 


template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesParameterJacobianEigen(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const scalar_t epsilon)
{
  auto numParams = ::pressio::ops::extent(mu,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  using dense_matrix_t = typename rom_data_t::dense_matrix_t;
  using vector_t = typename rom_data_t::vector_t;

  dense_matrix_t PhiTJPhi(romDim,numParams);
  PhiTJPhi.setZero();

  vector_t PhiTJPhiCol(romDim);
  PhiTJPhiCol.setZero();

  computeBasisTransposeTimesParameterJacobian_calc(appObj,y,mu,t,Phi, PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesParameterJacobianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const scalar_t epsilon)
{
  auto numParams = ::pressio::ops::extent(mu,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  using dense_matrix_t = typename rom_data_t::dense_matrix_t;
  using vector_t = typename rom_data_t::vector_t;

  dense_matrix_t PhiTJPhi("PhiJPphi_params",romDim,numParams);
  vector_t PhiTJPhiCol("PhiJPphiCol_params",romDim);

  computeBasisTransposeTimesParameterJacobian_calc(appObj,y,mu,t,Phi, PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 

//-----------------------------------



/*---------------------------------
Phi^T f computation
----------------------------------*/
template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
Eigen::Matrix<scalar_t,-1,1> computeBasisTransposeTimesVelocityEigen(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity(y,t,f);

  using vector_t = typename rom_data_t::vector_t;
  vector_t PhiTfCol(romDim);
  PhiTfCol.setZero();

  ::pressio::ops::product(::pressio::transpose(),1., Phi,f, 0.,PhiTfCol);
  return PhiTfCol;
} 

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesVelocityKokkos(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity(y,t,f);

  using vector_t		= typename rom_data_t::vector_t;
  vector_t PhiTfCol("PhiTfCol",romDim);

  ::pressio::ops::product(::pressio::transpose(),1., Phi,f, 0.,PhiTfCol);
  return PhiTfCol;
} 


// 
template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename dense_matrix_t, typename vector_t>
void computeBasisTransposeTimesJacobianTimesBasis_calc(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,dense_matrix_t & PhiTJPhi,vector_t & PhiTJPhiCol, const scalar_t epsilon){
  auto romDim = ::pressio::ops::extent(PhiTJPhi,0);

  auto ynorm = ::pressio::ops::norm2(y);
  // Create temporary working vectors
  auto fp = appObj.createVelocity();
  auto yp = ::pressio::ops::clone(y);
  auto workingVector = appObj.createVelocity();
  appObj.updateScalarParameters(mu);

  for  (int i = 0; i < romDim; i++){
    ::pressio::ops::set_zero(workingVector);

    auto PhiCol = getCol<scalar_t>(Phi,i);
    auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
    scalar_t eps = std::sqrt( (1. + ynorm)*epsilon ) / PhiColNorm;
    scalar_t scale = 0.5/eps;

    ::pressio::ops::update(yp,0.,y,1.,PhiCol,eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,0.,fp,1.*scale);

    ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,1.,fp,-1.*scale);

    ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., PhiTJPhiCol);
    setCol(PhiTJPhi,PhiTJPhiCol,i);
  }
}

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesJacobianTimesBasisKokkos(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  using vector_t		= typename rom_data_t::vector_t;
  using dense_matrix_t		= typename rom_data_t::dense_matrix_t;

  dense_matrix_t PhiTJPhi("PhiTJPhi",romDim,romDim); 
  vector_t PhiTJPhiCol("PhiTJPhiCol",romDim);
  computeBasisTransposeTimesJacobianTimesBasis_calc(appObj,y,mu,t,Phi,PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 


template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesJacobianTimesBasisEigen(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create result
  using vector_t		= typename rom_data_t::vector_t;
  using dense_matrix_t		= typename rom_data_t::dense_matrix_t;

  dense_matrix_t PhiTJPhi(romDim,romDim);
  PhiTJPhi.setZero();

  vector_t PhiTJPhiCol(romDim);
  PhiTJPhiCol.setZero();

  computeBasisTransposeTimesJacobianTimesBasis_calc(appObj,y,mu,t,Phi,PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeBasisTransposeTimesHessianTimesBasisTimesBasis_calc(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, hessian_t & BasisTransposeTimesHessianTimesBasisTimesBasis, vector_t & HCol,const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto ynorm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 
  for  (int i = 0; i < maxKForHessian; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        auto PhiCol = getCol<scalar_t>(Phi,i);
        auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
        scalar_t eps =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColNorm;
        auto scale = 1./(12.*eps*eps);
        //yp = y + 2.*eps*PhiCol;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol,2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.*scale);

        //yp = y + eps*PhiCol (or yp = yp - eps PhiCol);
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale);

    
        ::pressio::ops::update(workingVector,1.,f0,-30.*scale);
    
      
        //yp = y - eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale);
       
        //yp = y - 2.*eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
 
//        ::pressio::ops::scale(workingVector,1./(12.*(eps*eps)));
      }
      else{
        auto PhiCol_i = getCol<scalar_t>(Phi,i);
        auto PhiCol_j = getCol<scalar_t>(Phi,j);
        auto PhiColINorm = ::pressio::ops::norm2(PhiCol_i);
        auto PhiColJNorm = ::pressio::ops::norm2(PhiCol_j);
        scalar_t eps_i =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColINorm;
        scalar_t eps_j =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColJNorm;
        auto scale = 0.25/(eps_i*eps_j);


        //yp = y + eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i,PhiCol_j,eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);

        //yp = y + eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

        //yp = y - eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_i,-2.*eps_i,PhiCol_j,2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
       
        //yp = y - eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);
 
      }
      ::pressio::ops::product(::pressio::transpose(),1., Phi,workingVector , 0., HCol);
      setSymmetricHessianCol(BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,i,j);
    }
  }
} 



  
template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesHessianTimesBasisTimesBasisKokkos(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    =  typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesHessianTimesBasisTimesBasis("PhiTHPhiPhi", romDim,maxKForHessian,maxKForHessian);
  vector_t HCol("HCol",romDim);


  computeBasisTransposeTimesHessianTimesBasisTimesBasis_calc(appObj,y,mu,t,Phi,BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,maxKForHessian,epsilon);
  return BasisTransposeTimesHessianTimesBasisTimesBasis; 
} 



template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeBasisTransposeTimesHessianTimesBasisTimesBasisEigen(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const int maxKForHessian,const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    =  typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesHessianTimesBasisTimesBasis(romDim);

  vector_t HCol(romDim);
  HCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesHessianTimesBasisTimesBasis[i].resize(maxKForHessian);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < maxKForHessian; j++){
      BasisTransposeTimesHessianTimesBasisTimesBasis[i][j].resize(maxKForHessian);
      for (int k=0; k < maxKForHessian; k++){
        BasisTransposeTimesHessianTimesBasisTimesBasis[i][j][k] = 0.;
      }
    }
  }

  computeBasisTransposeTimesHessianTimesBasisTimesBasis_calc(appObj,y,mu,t,Phi,BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,maxKForHessian,epsilon);
  return BasisTransposeTimesHessianTimesBasisTimesBasis; 
} 

}}}
#endif
