#ifndef TPQL_MAIN_HPP_
#define TPQL_MAIN_HPP_
#include "tpql_rom_datatypes.hpp"
#include "tpql_operators.hpp"
#include "pressio/ops.hpp"

namespace pressio{ namespace rom{ namespace experimental{

template <typename rom_data_t, typename app_t, typename basis_t, typename scalar_t>
class TpwqRom
{
public:
  template <typename list_of_states_t, typename list_of_params_t, 
            typename list_of_times_t, typename list_of_coords_t, typename fom_state_t> 
  TpwqRom(const app_t & appObj,const basis_t & Phi,const list_of_states_t & InputListOfStatesAtLinearizationPoints,
          const list_of_params_t & InputListOfParamsAtLinearizationPoints, const list_of_times_t & InputListOfTimesAtLinearizationPoints, 
          const list_of_coords_t & InputListOfCoordsAtLinearizationPoints, const fom_state_t & fomReferenceState, const scalar_t & velocityErrorTolerance,
          const int maxKForHessian,
          const scalar_t epsilonForStateJacobian, const scalar_t epsilonForStateHessian, 
          const scalar_t epsilonForParameterJacobian, const scalar_t epsilonForParameterHessian, 
          const scalar_t epsilonForStateInMixedHessian, const scalar_t epsilonForParametersInMixedHessian) : Phi_(Phi) , appObj_(appObj), maxKForHessian_(maxKForHessian){

  romDim_ = ::pressio::ops::extent(Phi,1);
  auto param0 = InputListOfParamsAtLinearizationPoints[0];
  auto numParams = ::pressio::ops::extent(param0,0);

  // Compute mean of coordinates
  int numberOfInputPoints = InputListOfCoordsAtLinearizationPoints.size(); 
  coordinateDimension_ = InputListOfCoordsAtLinearizationPoints[0].size();

  std::vector<scalar_t> coordinateMeans(coordinateDimension_);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension_; j++){ 
      coordinateMeans[j] += InputListOfCoordsAtLinearizationPoints[i](j);
    }
  }

  for (int j=0; j < coordinateDimension_; j++){ 
   coordinateMeans[j] /= float(numberOfInputPoints);
  }

  std::vector<scalar_t> distanceFromMean(numberOfInputPoints);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension_; j++){ 
      distanceFromMean[i] += std::pow( InputListOfCoordsAtLinearizationPoints[i](j) - coordinateMeans[j] ,2.);
    }
    distanceFromMean[i] = std::sqrt(distanceFromMean[i]);
  }

  std::vector<int> sortedIndices(numberOfInputPoints);
  for (int i=0;i<numberOfInputPoints;i++){sortedIndices[i] = i;};
  std::sort( sortedIndices.begin(),sortedIndices.end(), [&](int i,int j){return distanceFromMean[i]<distanceFromMean[j];} );

  auto workingStateVector = appObj.createVelocity(); 

  int numCentroids = 0;
  for (int i = 0; i < numberOfInputPoints; i++){
    // Construct state at current point
    int sortedIndex = sortedIndices[i];
    auto currentState = InputListOfStatesAtLinearizationPoints[sortedIndex];
    auto currentTime = InputListOfTimesAtLinearizationPoints[sortedIndex];
    auto currentParam = InputListOfParamsAtLinearizationPoints[sortedIndex];
    auto currentCoords = InputListOfCoordsAtLinearizationPoints[sortedIndex];

    ::pressio::ops::product(::pressio::nontranspose(), 1., Phi,currentState,0.,workingStateVector);
    ::pressio::ops::update(workingStateVector,1.,fomReferenceState,1.); 
    if (i == 0){
      // Add point to coordinate list
      ListOfStatesAtLinearizationPoints_.push_back( currentState );  
      ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
      ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
      ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  


      auto PhiTf = ::pressio::rom::experimental::computeTestBasisTransposeTimesVelocityEigen<rom_data_t>(appObj, workingStateVector,currentParam,currentTime, Phi);
      ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);

      auto PhiTJPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesJacobianTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi,Phi,epsilonForStateJacobian);
      ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);

      auto PhiTHPhiPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesHessianTimesBasisTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi ,Phi, maxKForHessian,epsilonForStateHessian);
      ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);

      auto PhiTJParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterJacobianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , Phi, epsilonForParameterJacobian);
      ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);

      auto PhiTHParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi, Phi ,epsilonForParameterHessian);
      ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);

      auto PhiTHMixed = ::pressio::rom::experimental::computeTestBasisTransposeTimesMixedHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi,Phi, epsilonForStateHessian,epsilonForParametersInMixedHessian);
      ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);

      numCentroids += 1;      
    }
    else{

      auto relativeError = CheckGalerkinVelocity(currentState,workingStateVector,currentParam,currentTime,currentCoords);
      if (relativeError >= velocityErrorTolerance){
        auto fmt = ::pressio::utils::io::red() + ::pressio::utils::io::bold();
        ::pressio::utils::io::print_stdout(
          fmt, "Relative velocity error = ", relativeError , " adding new point " , ::pressio::utils::io::reset(), " \n");

        // Add point to coordinate list
        ListOfStatesAtLinearizationPoints_.push_back( currentState );  
        ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
        ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
        ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  
  

        auto PhiTf = ::pressio::rom::experimental::computeTestBasisTransposeTimesVelocityEigen<rom_data_t>(appObj, workingStateVector,currentParam,currentTime, Phi);
        ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);
  
        auto PhiTJPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesJacobianTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi ,Phi, epsilonForStateJacobian);
        ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);
  
        auto PhiTHPhiPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesHessianTimesBasisTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi,Phi, maxKForHessian,epsilonForStateHessian);
        ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);
  
        auto PhiTJParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterJacobianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi ,Phi, epsilonForParameterJacobian);
        ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);
  
        auto PhiTHParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , Phi, epsilonForParameterHessian);
        ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);
  
        auto PhiTHMixed = ::pressio::rom::experimental::computeTestBasisTransposeTimesMixedHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , Phi, epsilonForStateInMixedHessian,epsilonForParametersInMixedHessian);
        ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);
        numCentroids += 1;
      }
      else{
        auto fmt = ::pressio::utils::io::green() + ::pressio::utils::io::bold();
        ::pressio::utils::io::print_stdout(
          fmt, "Relative velocity error = ", relativeError , " skipping point " , ::pressio::utils::io::reset(), " \n");
      }
    }
  } 

  // Now export model to file
  std::ofstream outputFile ("tpwRomData.txt");
  if (outputFile.is_open())
  {
    // Write out basic information
    outputFile << romDim_ << std::endl;
    outputFile << maxKForHessian << std::endl;
    outputFile << numParams << std::endl;
    outputFile << coordinateDimension_ << std::endl;
    outputFile << numCentroids << std::endl;

    for (int c = 0; c < numCentroids; c++){

      // write out reduced states at centroids
      for (int i = 0; i < romDim_; i++){
        outputFile << ListOfStatesAtLinearizationPoints_[c](i) << std::endl;
      }
    }

      // write out parameters at centroids
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < numParams; i++){
        outputFile << ListOfParamsAtLinearizationPoints_[c](i) << std::endl;
      }
     }
      // write out times at centroids
    for (int c = 0; c < numCentroids; c++){
      outputFile << ListOfTimesAtLinearizationPoints_[c] << std::endl;
    }
      // write out coordinates at centroids
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < coordinateDimension_; i++){
        outputFile << ListOfCoordsAtLinearizationPoints_[c](i) << std::endl;
      }
     }
       // Write out PhiTf
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        outputFile << ListOfVelocitiesAtLinearizationPoints_[c](i) << std::endl;
      }
    }

       // Write out PhiTJPhi
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < romDim_; j++){
          outputFile << ListOfJacobiansAtLinearizationPoints_[c](i,j) << std::endl;
        }
      }
    }

      // Write out PhiTHPhiPhi
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < maxKForHessian; j++){
          for (int k = 0; k < maxKForHessian; k++){
            outputFile << ListOfHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }

      // Write out PhiTJParams
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < numParams; j++){
          outputFile << ListOfParamJacobiansAtLinearizationPoints_[c](i,j) << std::endl;
        }
      }
    }
      // Write out PhiTJParams
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < numParams; j++){
          for (int k = 0; k < numParams; k++){
            outputFile << ListOfParamHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }
      // Write out PhiTHMixed

    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < romDim_; j++){
          for (int k = 0; k < numParams; k++){
            outputFile << ListOfMixedHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }

    outputFile.close();
  }
  else std::cout << "Unable to open file" << std::endl;
 
  }
private:
  // types
  using rom_state_t = typename rom_data_t::vector_t;
  using params_t = typename rom_data_t::vector_t;
  using coords_t = typename rom_data_t::vector_t; ;

  using rom_velocity_t = typename rom_data_t::vector_t;
  using rom_jacobian_t = typename rom_data_t::dense_matrix_t;
  using rom_hessian_t = typename rom_data_t::dense_hessian_t; 

  using list_of_states_t = std::vector<rom_state_t>;
  using list_of_params_t = std::vector<params_t>;
  using list_of_times_t = std::vector<scalar_t>;
  using list_of_coords_t = std::vector<coords_t>;

  const app_t & appObj_;
  const basis_t & Phi_;
  int maxKForHessian_;
  int coordinateDimension_;
  int romDim_;
  list_of_states_t ListOfStatesAtLinearizationPoints_{0};
  list_of_params_t ListOfParamsAtLinearizationPoints_{0}; 
  list_of_times_t  ListOfTimesAtLinearizationPoints_{0};
  list_of_coords_t  ListOfCoordsAtLinearizationPoints_{0};

  std::vector< rom_velocity_t > ListOfVelocitiesAtLinearizationPoints_{0};
  std::vector< rom_jacobian_t > ListOfJacobiansAtLinearizationPoints_{0};
  std::vector< rom_hessian_t  > ListOfHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_jacobian_t  > ListOfParamJacobiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfParamHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfMixedHessiansAtLinearizationPoints_{0}; 

public:
  template <typename params_t, typename coords_t>
  rom_state_t GalerkinVelocity(const rom_state_t & xhat,const params_t & mu,const scalar_t t,const coords_t & PointCoords){
    int linPointIndx = FindClosestLinearizationPoint(PointCoords);
    auto VelocityAtLinearizationPoint = ListOfVelocitiesAtLinearizationPoints_[linPointIndx];
    auto JacobianAtLinearizationPoint = ListOfJacobiansAtLinearizationPoints_[linPointIndx];
    auto HessianAtLinearizationPoint =  ListOfHessiansAtLinearizationPoints_[linPointIndx];
    auto ParameterJacobianAtLinearizationPoint = ListOfParamJacobiansAtLinearizationPoints_[linPointIndx];
    auto ParameterHessianAtLinearizationPoint = ListOfParamHessiansAtLinearizationPoints_[linPointIndx];
    auto MixedHessianAtLinearizationPoint = ListOfMixedHessiansAtLinearizationPoints_[linPointIndx];
    auto ReducedStateAtLinearizationPoint = ListOfStatesAtLinearizationPoints_[linPointIndx];
    auto ParamAtLinearizationPoint = ListOfParamsAtLinearizationPoints_[linPointIndx];
    auto CoordsAtLinearizationPoint = ListOfCoordsAtLinearizationPoints_[linPointIndx];
    auto fmt1 = ::pressio::utils::io::cyan() + ::pressio::utils::io::bold();

    ::pressio::utils::io::print_stdout(
        fmt1, "Calling Galerkin velocity with coordinates = ");
    auto numCoords = ::pressio::ops::extent(PointCoords,0);
    for (int i = 0; i < numCoords; i++){
      ::pressio::utils::io::print_stdout(fmt1, " " , PointCoords(i) );
    }
    ::pressio::utils::io::print_stdout("\n");
    ::pressio::utils::io::print_stdout(fmt1, "Closest linearization point = ");
    for (int i = 0; i < numCoords; i++){
      ::pressio::utils::io::print_stdout(fmt1, " " , CoordsAtLinearizationPoint(i) );
    }
    ::pressio::utils::io::print_stdout("\n");


    rom_state_t d_xHat(xhat);
    d_xHat.setZero();

    rom_state_t d_xHat_Reduced(maxKForHessian_);
    d_xHat_Reduced.setZero();

    ::pressio::ops::update(d_xHat,0.,xhat,1.,ReducedStateAtLinearizationPoint,-1.);

    for (int i = 0; i < maxKForHessian_; i++){
      d_xHat_Reduced(i) = d_xHat(i);
    }

    params_t d_mu(mu);
    ::pressio::ops::set_zero(d_mu);
    ::pressio::ops::update(d_mu,0.,mu,1.,ParamAtLinearizationPoint,-1.);

    rom_velocity_t f(romDim_);
    ::pressio::ops::set_zero(f);
    auto budget0 = ::pressio::ops::norm2(f);
    ::pressio::ops::update(f,0.,VelocityAtLinearizationPoint,1.);
    auto budget1 = ::pressio::ops::norm2(f);
    auto budget1Diff = std::sqrt( std::abs(budget1*budget1 - budget0*budget0) + 1.e-30);

    ::pressio::ops::product(::pressio::nontranspose(),1., JacobianAtLinearizationPoint,         d_xHat,                    1.,f);

    auto budget2 = ::pressio::ops::norm2(f);
    auto budget2Diff = std::sqrt( std::abs(budget2*budget2 - budget1*budget1) + 1.e-30);

    ::pressio::ops::product(::pressio::nontranspose(),1., ParameterJacobianAtLinearizationPoint,d_mu  ,                    1.,f);

    auto budget3 = ::pressio::ops::norm2(f);
    auto budget3Diff = std::sqrt( std::abs(budget3*budget3 - budget2*budget2) + 1.e-30);
    
    tensorMultiply(HessianAtLinearizationPoint, d_xHat_Reduced,d_xHat_Reduced,0.5,f); 

    auto budget4 = ::pressio::ops::norm2(f);
    auto budget4Diff = std::sqrt( std::abs(budget4*budget4 - budget3*budget3) + 1.e-30);

    tensorMultiply(ParameterHessianAtLinearizationPoint, d_mu,d_mu,0.5,f); 

    auto budget5 = ::pressio::ops::norm2(f);
    auto budget5Diff = std::sqrt( std::abs(budget5*budget5 - budget4*budget4) + 1.e-30);

    tensorMultiply(MixedHessianAtLinearizationPoint, d_xHat,d_mu,1.,f); 

    auto budget6 = ::pressio::ops::norm2(f);
    auto budget6Diff = std::sqrt( std::abs(budget6*budget6 - budget5*budget5) + 1.e-30);

    ::pressio::utils::io::print_stdout("\n");
    auto fmt = ::pressio::utils::io::blue() + ::pressio::utils::io::bold();
    ::pressio::utils::io::print_stdout(
        fmt, "Budget 0 = ", budget0 , ::pressio::utils::io::reset(), " | ");
    ::pressio::utils::io::print_stdout(
        fmt, "linearization point budget = ", budget1Diff , ::pressio::utils::io::reset(), " | ");
    ::pressio::utils::io::print_stdout(
        fmt, "State Jacobian budget = ", budget2Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Parameter Jacobian budget = ", budget3Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "State Hessian budget = ", budget4Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Parameter Hessian budget = ", budget5Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Mixed Hessian budget = ", budget6Diff , ::pressio::utils::io::reset(), "\n");

    return f;
  }

private:
  void tensorMultiply(const rom_hessian_t & H,const rom_state_t & a1,const rom_state_t & a2,const scalar_t alpha, rom_velocity_t & f){
    auto d1 = H.size();
    auto d2 = H[0].size();
    auto d3 = H[0][0].size();
    for (int i = 0; i < d1; i++){
      for (int j = 0; j < d2; j++){
        for (int k = 0; k < d3; k++){
          f[i] += alpha*H[i][j][k] * a1(j) * a2(k);
        }
      }
    } 
  }

  template <typename fom_state_t, typename params_t, typename coords_t>
  double CheckGalerkinVelocity(const rom_state_t & xhat,const fom_state_t & x, const params_t & mu,const scalar_t t,const coords_t & PointCoords)
  {
    auto fTpwq = GalerkinVelocity(xhat,mu,t,PointCoords);
    auto fTrue = appObj_.createVelocity();
    ::pressio::ops::set_zero(fTrue);
    appObj_.updateScalarParameters(mu);
    appObj_.velocity(x,t,fTrue);


    rom_velocity_t fGalerkin(romDim_);
    ::pressio::ops::set_zero(fGalerkin);
    rom_velocity_t velocityError(romDim_);
    ::pressio::ops::set_zero(velocityError);
    ::pressio::ops::product(::pressio::transpose(),1., Phi_,fTrue, 0.,fGalerkin);    
    ::pressio::ops::update(velocityError,0.,fTpwq,1.,fGalerkin,-1.);


    auto galerkinNorm = ::pressio::ops::norm2(fGalerkin);
    auto tpwNorm = ::pressio::ops::norm2(fTpwq);
    auto fmt = ::pressio::utils::io::blue() + ::pressio::utils::io::bold();
    ::pressio::utils::io::print_stdout(
        fmt, "Galerkin velocity norm = ", galerkinNorm , " tpw velocity norm = " , tpwNorm , ::pressio::utils::io::reset(), "\n");

    auto errorNorm = ::pressio::ops::norm2(velocityError)/ ::pressio::ops::norm2(fGalerkin); 
    return errorNorm; 
  }

  template<typename coords_t>
  int FindClosestLinearizationPoint(coords_t PointCoords){

    int numberOfAddedLinearizationPoints = ListOfCoordsAtLinearizationPoints_.size();
    std::vector<scalar_t> distance(numberOfAddedLinearizationPoints);

    for (int i = 0; i < numberOfAddedLinearizationPoints; i++){
      for (int j=0; j < coordinateDimension_; j++){
        distance[i] += std::pow( ListOfCoordsAtLinearizationPoints_[i](j) - PointCoords[j] ,2.);
      }
      distance[i] = std::sqrt(distance[i]);
    }
    int linPointIndex = std::min_element(distance.begin(), distance.end()) - distance.begin();
    return linPointIndex;
  }

  //template <typename params_t, typename coords_t>
  //void GalerkinVelocity(const rom_state_t & xhat,const params_t & mu,const scalar_t t,const coords_t & PointCoords){
  //}
};


template <typename rom_data_t, typename app_t, typename basis_t, typename scalar_t>
class TpwqQuasiLspgRom
{
public:
  template <typename list_of_states_t, typename list_of_params_t, 
            typename list_of_times_t, typename list_of_coords_t, typename fom_state_t> 
  TpwqQuasiLspgRom(const app_t & appObj,const basis_t & Phi,const list_of_states_t & InputListOfStatesAtLinearizationPoints,
          const list_of_params_t & InputListOfParamsAtLinearizationPoints, const list_of_times_t & InputListOfTimesAtLinearizationPoints, 
          const list_of_coords_t & InputListOfCoordsAtLinearizationPoints, const fom_state_t & fomReferenceState, const scalar_t & velocityErrorTolerance,
          const int maxKForHessian,
          const scalar_t epsilonForStateJacobian, const scalar_t epsilonForStateHessian, 
          const scalar_t epsilonForParameterJacobian, const scalar_t epsilonForParameterHessian, 
          const scalar_t epsilonForStateInMixedHessian, const scalar_t epsilonForParametersInMixedHessian) : appObj_(appObj), Phi_(Phi),  maxKForHessian_(maxKForHessian){

  romDim_ = ::pressio::ops::extent(Phi,1);
  auto param0 = InputListOfParamsAtLinearizationPoints[0];
  auto numParams = ::pressio::ops::extent(param0,0);

  // Compute mean of coordinates
  int numberOfInputPoints = InputListOfCoordsAtLinearizationPoints.size(); 
  coordinateDimension_ = InputListOfCoordsAtLinearizationPoints[0].size();

  std::vector<scalar_t> coordinateMeans(coordinateDimension_);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension_; j++){ 
      coordinateMeans[j] += InputListOfCoordsAtLinearizationPoints[i](j);
    }
  }

  for (int j=0; j < coordinateDimension_; j++){ 
   coordinateMeans[j] /= float(numberOfInputPoints);
  }

  std::vector<scalar_t> distanceFromMean(numberOfInputPoints);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension_; j++){ 
      distanceFromMean[i] += std::pow( InputListOfCoordsAtLinearizationPoints[i](j) - coordinateMeans[j] ,2.);
    }
    distanceFromMean[i] = std::sqrt(distanceFromMean[i]);
  }

  std::vector<int> sortedIndices(numberOfInputPoints);
  for (int i=0;i<numberOfInputPoints;i++){sortedIndices[i] = i;};
  std::sort( sortedIndices.begin(),sortedIndices.end(), [&](int i,int j){return distanceFromMean[i]<distanceFromMean[j];} );

  auto workingStateVector = appObj.createVelocity(); 

  int numCentroids = 0;
  for (int i = 0; i < numberOfInputPoints; i++){
    // Construct state at current point
    int sortedIndex = sortedIndices[i];
    auto currentState = InputListOfStatesAtLinearizationPoints[sortedIndex];
    auto currentTime = InputListOfTimesAtLinearizationPoints[sortedIndex];
    auto currentParam = InputListOfParamsAtLinearizationPoints[sortedIndex];
    auto currentCoords = InputListOfCoordsAtLinearizationPoints[sortedIndex];

    ::pressio::ops::product(::pressio::nontranspose(), 1., Phi,currentState,0.,workingStateVector);
    ::pressio::ops::update(workingStateVector,1.,fomReferenceState,1.); 
    if (i == 0){
      // Add point to coordinate list
      ListOfStatesAtLinearizationPoints_.push_back( currentState );  
      ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
      ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
      ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  

      //std::cout << " 1 " <<std::endl;
      auto JPhi = ::pressio::ops::clone(Phi);
      ::pressio::rom::experimental::computeJacobianTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi,JPhi,epsilonForStateJacobian);
      ListOfJPhisAtLinearizationPoints_.push_back(JPhi);
      //std::cout << " 2 " <<std::endl;

      auto PhiTf = ::pressio::rom::experimental::computeTestBasisTransposeTimesVelocityEigen<rom_data_t>(appObj, workingStateVector,currentParam,currentTime,  JPhi);
      ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);
      //std::cout << " 3 " <<std::endl;

      auto JPhiVecView2 = JPhi.getMultiVectorView();

      //std::cout << "JPhiVewView in main " << JPhiVecView2.getLocalLength() << std::endl;

      auto PhiTJPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesJacobianTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi,JPhi,epsilonForStateJacobian);
      ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);
      //std::cout << " 4 " <<std::endl;

      auto PhiTHPhiPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesHessianTimesBasisTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , JPhi, maxKForHessian,epsilonForStateHessian);
      ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);
      //std::cout << " 5 " <<std::endl;

      auto PhiTJParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterJacobianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi ,JPhi,epsilonForParameterJacobian);
      ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);
      //std::cout << " 6 " <<std::endl;

      auto PhiTHParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi, JPhi, epsilonForParameterHessian);
      ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);
      //std::cout << " 7 " <<std::endl;

      auto PhiTHMixed = ::pressio::rom::experimental::computeTestBasisTransposeTimesMixedHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi,JPhi, epsilonForStateHessian,epsilonForParametersInMixedHessian);
      ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);
      //std::cout << " 8 " <<std::endl;

      numCentroids += 1;      
    }
    else{

      auto relativeError = CheckQuasiLspgVelocity(currentState,workingStateVector,currentParam,currentTime,currentCoords);
      if (relativeError >= velocityErrorTolerance){
        auto fmt = ::pressio::utils::io::red() + ::pressio::utils::io::bold();
        ::pressio::utils::io::print_stdout(
          fmt, "Relative velocity error = ", relativeError , " adding new point " , ::pressio::utils::io::reset(), " \n");

        // Add point to coordinate list
        ListOfStatesAtLinearizationPoints_.push_back( currentState );  
        ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
        ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
        ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  
  
        auto JPhi = ::pressio::ops::clone(Phi);
        ::pressio::rom::experimental::computeJacobianTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi,JPhi,epsilonForStateJacobian);
        ListOfJPhisAtLinearizationPoints_.push_back(JPhi);

        auto PhiTf = ::pressio::rom::experimental::computeTestBasisTransposeTimesVelocityEigen<rom_data_t>(appObj, workingStateVector,currentParam,currentTime, JPhi);
        ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);
  
        auto PhiTJPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesJacobianTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , JPhi, epsilonForStateJacobian);
        ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);
  
        auto PhiTHPhiPhi = ::pressio::rom::experimental::computeTestBasisTransposeTimesHessianTimesBasisTimesBasisEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi, JPhi, maxKForHessian,epsilonForStateHessian);
        ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);
  
        auto PhiTJParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterJacobianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , JPhi, epsilonForParameterJacobian);
        ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);
  
        auto PhiTHParams = ::pressio::rom::experimental::computeTestBasisTransposeTimesParameterHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi , JPhi, epsilonForParameterHessian);
        ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);
  
        auto PhiTHMixed = ::pressio::rom::experimental::computeTestBasisTransposeTimesMixedHessianEigen<rom_data_t>(appObj,workingStateVector,currentParam,currentTime, Phi ,JPhi,epsilonForStateInMixedHessian,epsilonForParametersInMixedHessian);
        ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);
        numCentroids += 1;
      }
      else{
        auto fmt = ::pressio::utils::io::green() + ::pressio::utils::io::bold();
        ::pressio::utils::io::print_stdout(
          fmt, "Relative velocity error = ", relativeError , " skipping point " , ::pressio::utils::io::reset(), " \n");
      }
    }
  } 

  // Now export model to file
  std::ofstream outputFile ("tpwRomData.txt");
  if (outputFile.is_open())
  {
    // Write out basic information
    outputFile << romDim_ << std::endl;
    outputFile << maxKForHessian << std::endl;
    outputFile << numParams << std::endl;
    outputFile << coordinateDimension_ << std::endl;
    outputFile << numCentroids << std::endl;

    for (int c = 0; c < numCentroids; c++){

      // write out reduced states at centroids
      for (int i = 0; i < romDim_; i++){
        outputFile << ListOfStatesAtLinearizationPoints_[c](i) << std::endl;
      }
    }

      // write out parameters at centroids
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < numParams; i++){
        outputFile << ListOfParamsAtLinearizationPoints_[c](i) << std::endl;
      }
     }
      // write out times at centroids
    for (int c = 0; c < numCentroids; c++){
      outputFile << ListOfTimesAtLinearizationPoints_[c] << std::endl;
    }
      // write out coordinates at centroids
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < coordinateDimension_; i++){
        outputFile << ListOfCoordsAtLinearizationPoints_[c](i) << std::endl;
      }
     }
       // Write out PhiTf
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        outputFile << ListOfVelocitiesAtLinearizationPoints_[c](i) << std::endl;
      }
    }

       // Write out PhiTJPhi
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < romDim_; j++){
          outputFile << ListOfJacobiansAtLinearizationPoints_[c](i,j) << std::endl;
        }
      }
    }

      // Write out PhiTHPhiPhi
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < maxKForHessian; j++){
          for (int k = 0; k < maxKForHessian; k++){
            outputFile << ListOfHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }

      // Write out PhiTJParams
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < numParams; j++){
          outputFile << ListOfParamJacobiansAtLinearizationPoints_[c](i,j) << std::endl;
        }
      }
    }
      // Write out PhiTJParams
    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < numParams; j++){
          for (int k = 0; k < numParams; k++){
            outputFile << ListOfParamHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }
      // Write out PhiTHMixed

    for (int c = 0; c < numCentroids; c++){
      for (int i = 0; i < romDim_; i++){
        for (int j = 0; j < romDim_; j++){
          for (int k = 0; k < numParams; k++){
            outputFile << ListOfMixedHessiansAtLinearizationPoints_[c][i][j][k] << std::endl;
          }
        }
      }
    }

    outputFile.close();
  }
  else std::cout << "Unable to open file" << std::endl;
 
  }
private:
  // types
  using rom_state_t = typename rom_data_t::vector_t;
  using params_t = typename rom_data_t::vector_t;
  using coords_t = typename rom_data_t::vector_t; ;

  using rom_velocity_t = typename rom_data_t::vector_t;
  using rom_jacobian_t = typename rom_data_t::dense_matrix_t;
  using rom_hessian_t = typename rom_data_t::dense_hessian_t; 

  using list_of_states_t = std::vector<rom_state_t>;
  using list_of_params_t = std::vector<params_t>;
  using list_of_times_t = std::vector<scalar_t>;
  using list_of_coords_t = std::vector<coords_t>;
  using list_of_jphi_t = std::vector<basis_t>;

  const app_t & appObj_;
  const basis_t & Phi_;
  int maxKForHessian_;
  int coordinateDimension_;
  int romDim_;
  list_of_states_t ListOfStatesAtLinearizationPoints_{0};
  list_of_params_t ListOfParamsAtLinearizationPoints_{0}; 
  list_of_times_t  ListOfTimesAtLinearizationPoints_{0};
  list_of_coords_t  ListOfCoordsAtLinearizationPoints_{0};
  list_of_jphi_t ListOfJPhisAtLinearizationPoints_{0}; 

  std::vector< rom_velocity_t > ListOfVelocitiesAtLinearizationPoints_{0};
  std::vector< rom_jacobian_t > ListOfJacobiansAtLinearizationPoints_{0};
  std::vector< rom_hessian_t  > ListOfHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_jacobian_t  > ListOfParamJacobiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfParamHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfMixedHessiansAtLinearizationPoints_{0}; 

public:
  template <typename params_t, typename coords_t>
  rom_state_t QuasiLspgVelocity(const rom_state_t & xhat,const params_t & mu,const scalar_t t,const coords_t & PointCoords){
    int linPointIndx = FindClosestLinearizationPoint(PointCoords);
    auto VelocityAtLinearizationPoint = ListOfVelocitiesAtLinearizationPoints_[linPointIndx];
    auto JacobianAtLinearizationPoint = ListOfJacobiansAtLinearizationPoints_[linPointIndx];
    auto HessianAtLinearizationPoint =  ListOfHessiansAtLinearizationPoints_[linPointIndx];
    auto ParameterJacobianAtLinearizationPoint = ListOfParamJacobiansAtLinearizationPoints_[linPointIndx];
    auto ParameterHessianAtLinearizationPoint = ListOfParamHessiansAtLinearizationPoints_[linPointIndx];
    auto MixedHessianAtLinearizationPoint = ListOfMixedHessiansAtLinearizationPoints_[linPointIndx];
    auto ReducedStateAtLinearizationPoint = ListOfStatesAtLinearizationPoints_[linPointIndx];
    auto ParamAtLinearizationPoint = ListOfParamsAtLinearizationPoints_[linPointIndx];
    auto CoordsAtLinearizationPoint = ListOfCoordsAtLinearizationPoints_[linPointIndx];
    auto fmt1 = ::pressio::utils::io::cyan() + ::pressio::utils::io::bold();

    ::pressio::utils::io::print_stdout(
        fmt1, "Calling quasi LSPG velocity with coordinates = ");
    auto numCoords = ::pressio::ops::extent(PointCoords,0);
    for (int i = 0; i < numCoords; i++){
      ::pressio::utils::io::print_stdout(fmt1, " " , PointCoords(i) );
    }
    ::pressio::utils::io::print_stdout("\n");
    ::pressio::utils::io::print_stdout(fmt1, "Closest linearization point = ");
    for (int i = 0; i < numCoords; i++){
      ::pressio::utils::io::print_stdout(fmt1, " " , CoordsAtLinearizationPoint(i) );
    }
    ::pressio::utils::io::print_stdout("\n");


    rom_state_t d_xHat(xhat);
    d_xHat.setZero();

    rom_state_t d_xHat_Reduced(maxKForHessian_);
    d_xHat_Reduced.setZero();

    ::pressio::ops::update(d_xHat,0.,xhat,1.,ReducedStateAtLinearizationPoint,-1.);

    for (int i = 0; i < maxKForHessian_; i++){
      d_xHat_Reduced(i) = d_xHat(i);
    }

    params_t d_mu(mu);
    ::pressio::ops::set_zero(d_mu);
    ::pressio::ops::update(d_mu,0.,mu,1.,ParamAtLinearizationPoint,-1.);

    rom_velocity_t f(romDim_);
    ::pressio::ops::set_zero(f);
    auto budget0 = ::pressio::ops::norm2(f);
    ::pressio::ops::update(f,0.,VelocityAtLinearizationPoint,1.);
    auto budget1 = ::pressio::ops::norm2(f);
    auto budget1Diff = std::sqrt( std::abs(budget1*budget1 - budget0*budget0) + 1.e-30);

    ::pressio::ops::product(::pressio::nontranspose(),1., JacobianAtLinearizationPoint,         d_xHat,                    1.,f);

    auto budget2 = ::pressio::ops::norm2(f);
    auto budget2Diff = std::sqrt( std::abs(budget2*budget2 - budget1*budget1) + 1.e-30);

    ::pressio::ops::product(::pressio::nontranspose(),1., ParameterJacobianAtLinearizationPoint,d_mu  ,                    1.,f);

    auto budget3 = ::pressio::ops::norm2(f);
    auto budget3Diff = std::sqrt( std::abs(budget3*budget3 - budget2*budget2) + 1.e-30);
    
    tensorMultiply(HessianAtLinearizationPoint, d_xHat_Reduced,d_xHat_Reduced,0.5,f); 

    auto budget4 = ::pressio::ops::norm2(f);
    auto budget4Diff = std::sqrt( std::abs(budget4*budget4 - budget3*budget3) + 1.e-30);

    tensorMultiply(ParameterHessianAtLinearizationPoint, d_mu,d_mu,0.5,f); 

    auto budget5 = ::pressio::ops::norm2(f);
    auto budget5Diff = std::sqrt( std::abs(budget5*budget5 - budget4*budget4) + 1.e-30);

    tensorMultiply(MixedHessianAtLinearizationPoint, d_xHat,d_mu,1.,f); 

    auto budget6 = ::pressio::ops::norm2(f);
    auto budget6Diff = std::sqrt( std::abs(budget6*budget6 - budget5*budget5) + 1.e-30);

    ::pressio::utils::io::print_stdout("\n");
    auto fmt = ::pressio::utils::io::blue() + ::pressio::utils::io::bold();
    ::pressio::utils::io::print_stdout(
        fmt, "Budget 0 = ", budget0 , ::pressio::utils::io::reset(), " | ");
    ::pressio::utils::io::print_stdout(
        fmt, "linearization point budget = ", budget1Diff , ::pressio::utils::io::reset(), " | ");
    ::pressio::utils::io::print_stdout(
        fmt, "State Jacobian budget = ", budget2Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Parameter Jacobian budget = ", budget3Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "State Hessian budget = ", budget4Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Parameter Hessian budget = ", budget5Diff , ::pressio::utils::io::reset(), " | ");

    ::pressio::utils::io::print_stdout(
        fmt, "Mixed Hessian budget = ", budget6Diff , ::pressio::utils::io::reset(), "\n");

    return f;
  }

private:
  void tensorMultiply(const rom_hessian_t & H,const rom_state_t & a1,const rom_state_t & a2,const scalar_t alpha, rom_velocity_t & f){
    auto d1 = H.size();
    auto d2 = H[0].size();
    auto d3 = H[0][0].size();
    for (int i = 0; i < d1; i++){
      for (int j = 0; j < d2; j++){
        for (int k = 0; k < d3; k++){
          f[i] += alpha*H[i][j][k] * a1(j) * a2(k);
        }
      }
    } 
  }

  template <typename fom_state_t, typename params_t, typename coords_t>
  double CheckQuasiLspgVelocity(const rom_state_t & xhat,const fom_state_t & x, const params_t & mu,const scalar_t t,const coords_t & PointCoords)
  {
    auto fTpwq = QuasiLspgVelocity(xhat,mu,t,PointCoords);
    auto fTrue = appObj_.createVelocity();
    ::pressio::ops::set_zero(fTrue);
    appObj_.updateScalarParameters(mu);
    appObj_.velocity(x,t,fTrue);

    int linPointIndx = FindClosestLinearizationPoint(PointCoords);
    auto JPhiAtLinearizationPoint = ListOfJPhisAtLinearizationPoints_[linPointIndx];

    rom_velocity_t fQuasiLspg(romDim_);
    ::pressio::ops::set_zero(fQuasiLspg);
    rom_velocity_t velocityError(romDim_);
    ::pressio::ops::set_zero(velocityError);
    ::pressio::ops::product(::pressio::transpose(),1., JPhiAtLinearizationPoint,fTrue, 0.,fQuasiLspg);    
    ::pressio::ops::update(velocityError,0.,fTpwq,1.,fQuasiLspg,-1.);


    auto QuasiLspgNorm = ::pressio::ops::norm2(fQuasiLspg);
    auto tpwNorm = ::pressio::ops::norm2(fTpwq);
    auto fmt = ::pressio::utils::io::blue() + ::pressio::utils::io::bold();
    ::pressio::utils::io::print_stdout(
        fmt, "Quasi LSPG velocity norm = ", QuasiLspgNorm , " tpw velocity norm = " , tpwNorm , ::pressio::utils::io::reset(), "\n");

    auto errorNorm = ::pressio::ops::norm2(velocityError)/ ::pressio::ops::norm2(fQuasiLspg); 
    return errorNorm; 
  }

  template<typename coords_t>
  int FindClosestLinearizationPoint(coords_t PointCoords){

    int numberOfAddedLinearizationPoints = ListOfCoordsAtLinearizationPoints_.size();
    std::vector<scalar_t> distance(numberOfAddedLinearizationPoints);

    for (int i = 0; i < numberOfAddedLinearizationPoints; i++){
      for (int j=0; j < coordinateDimension_; j++){
        distance[i] += std::pow( ListOfCoordsAtLinearizationPoints_[i](j) - PointCoords[j] ,2.);
      }
      distance[i] = std::sqrt(distance[i]);
    }
    int linPointIndex = std::min_element(distance.begin(), distance.end()) - distance.begin();
    return linPointIndex;
  }

};


}}}
#endif
