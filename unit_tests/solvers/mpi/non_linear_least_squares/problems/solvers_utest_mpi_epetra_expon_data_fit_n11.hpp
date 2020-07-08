
#ifndef SOLVERS_UTEST_MPI_EPETRA_EXPONENTIAL_DATA_FIT_N11_HPP_
#define SOLVERS_UTEST_MPI_EPETRA_EXPONENTIAL_DATA_FIT_N11_HPP_

#include "Epetra_MpiComm.h"
#include "pressio_solvers.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.5
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f
 */

namespace pressio{ namespace solvers{ namespace test{

struct ExpDataFitN11
{
  const int numEq_ = 65;
  const int numUn_ = 11;

  template <typename T>
  using shptr = std::shared_ptr<T>;

  using state_type    = pressio::containers::Vector<Eigen::VectorXd>;
  using residual_type = pressio::containers::Vector<Epetra_Vector>;
  using nat_vec_type  = pressio::containers::Vector<Epetra_Vector>;
  using jacobian_type = pressio::containers::MultiVector<Epetra_MultiVector>;
  using scalar_type   = double;

  shptr<Epetra_MpiComm> comm_ = {};
  int rank_		      = {};
  int numProc_		      = {};

  shptr<Epetra_Map> rowMap_   = {};
  int NumMyElem_	      = {}; // num of my elements
  std::vector<int> myGel_     = {}; // my global elements

  shptr<residual_type> R_     = {};
  shptr<jacobian_type> J_     = {};
  shptr<nat_vec_type> tt_     = {}; // store the t value for the data
  shptr<nat_vec_type> yy_     = {}; // store the y values

  const std::vector<double> trueS = {1.3099771546537, 0.43155379466611,
				     0.63366169898352, 0.59943053483825,
				     0.75418322650124, 0.90428857933176,
				     1.3658118355256,  4.82369881684,
				     2.3986848661039,  4.5688745976743,
				     5.67534147057};

  void storeTimes(){
    int i=0;
    for (auto const & it : myGel_){
      (*tt_)[i] = static_cast<scalar_type>(it)/10.;
      i++;
    };
  }

  void storeYValues(){
    const std::vector<double> allYs = {1.366, 1.191, 1.112, 1.013,
				       9.91e-1, 8.85e-1, 8.31e-1, 8.47e-1,
				       7.86e-1, 7.25e-1, 7.46e-1, 6.79e-1,
				       6.08e-1, 6.55e-1, 6.16e-1, 6.06e-1,
				       6.02e-1, 6.26e-1, 6.51e-1, 7.24e-1,
				       6.49e-1, 6.49e-1, 6.94e-1, 6.44e-1,
				       6.24e-1, 6.61e-1, 6.12e-1, 5.58e-1,
				       5.33e-1, 4.95e-1, 5.0e-1,  4.23e-1,
				       3.95e-1, 3.75e-1, 3.72e-1, 3.91e-1,
				       3.96e-1, 4.05e-1, 4.28e-1, 4.29e-1,
				       5.23e-1, 5.62e-1, 6.07e-1, 6.53e-1,
				       6.72e-1, 7.08e-1, 6.33e-1, 6.68e-1,
				       6.45e-1, 6.32e-1, 5.91e-1, 5.59e-1,
				       5.97e-1, 6.25e-1, 7.39e-1, 7.1e-1,
				       7.29e-1, 7.2e-1, 6.36e-1, 5.81e-1,
				       4.28e-1, 2.92e-1, 1.62e-1, 9.8e-2, 5.4e-2};

    int i=0;
    for (auto const & it : myGel_){
      (*yy_)[i] = allYs[it];
      i++;
    };
  }

  ExpDataFitN11(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = std::make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    assert(numProc_==2 or numProc_==3);

    // create map
    rowMap_ = std::make_shared<Epetra_Map>(numEq_, 0, *comm_);
    // store my elements
    NumMyElem_ = rowMap_->NumMyElements();
    myGel_.resize(NumMyElem_);
    rowMap_->MyGlobalElements(myGel_.data());

    R_ = std::make_shared<residual_type>(*rowMap_);
    J_ = std::make_shared<jacobian_type>(*rowMap_, numUn_);
    tt_ = std::make_shared<nat_vec_type>(*rowMap_);
    yy_ = std::make_shared<nat_vec_type>(*rowMap_);
    storeTimes();
    storeYValues();
    std::cout << std::endl;
  }//setUp

  inline scalar_type model(const state_type & x, scalar_type t)const{
    auto temp1 = exp(-x[4]*t);
    auto temp2 = exp(-x[5] * (t-x[8]) * (t-x[8]) );
    auto temp3 = exp(-x[6] * (t-x[9]) * (t-x[9]) );
    auto temp4 = exp(-x[7] * (t-x[10]) * (t-x[10]) );
    return x[0]*temp1 + x[1]*temp2 + x[2]*temp3 + x[3]*temp4;
  }

  residual_type createResidualObject(const state_type& x) const {
    // residual(x, *R_);
    return *R_;
  }

  jacobian_type createJacobianObject(const state_type & x) const{
    // jacobian(x, *J_);
    return *J_;
  }//end jacobian

  void residual(const state_type& x, residual_type & R,
    ::pressio::Norm normKind,
    scalar_type & normResidual) const 
  {
    for (auto i=0; i< NumMyElem_; i++){
      R[i] = (*yy_)[i] - this->model(x, (*tt_)[i]);
    };

    if (normKind == pressio::Norm::L2) R.data()->Norm2(&normResidual);
    if (normKind == pressio::Norm::L1) R.data()->Norm1(&normResidual);
  }

  void jacobian(const state_type & x, jacobian_type & jac) const{
    for (int i=0; i<NumMyElem_; i++)
      {
	scalar_type t = (*tt_)[i];

	auto temp1 = exp(-x[4]*t);
	auto temp2 = exp(-x[5] * (t-x[8]) * (t-x[8]) );
	auto temp3 = exp(-x[6] * (t-x[9]) * (t-x[9]) );
	auto temp4 = exp(-x[7] * (t-x[10]) * (t-x[10]) );
	jac(i,0) = -temp1;
	jac(i,1) = -temp2;
	jac(i,2) = -temp3;
	jac(i,3) = -temp4;
	jac(i,4) = x[0]*temp1*t;
	jac(i,5) = x[1] * (t-x[8]) * (t-x[8])  * temp2;
	jac(i,6) = x[2] * (t-x[9]) * (t-x[9])  * temp3;
	jac(i,7) = x[3] * (t-x[10]) * (t-x[10]) * temp4;
	jac(i,8) =  -2.*x[1]*x[5]*(t-x[8])*temp2;
	jac(i,9) =  -2.*x[2]*x[6]*(t-x[9])*temp3;
	jac(i,10) = -2.*x[3]*x[7]*(t-x[10])*temp4;
      }
    //    jac.data()->Print(std::cout);
  }//end jacobian


};

}}} //end namespace pressio::solvers::test
#endif
