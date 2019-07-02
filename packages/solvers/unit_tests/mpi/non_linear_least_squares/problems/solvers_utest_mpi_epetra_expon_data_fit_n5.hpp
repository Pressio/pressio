
#ifndef SOLVERS_UTEST_MPI_EPETRA_EXPONENTIAL_DATA_FIT_N5_HPP_
#define SOLVERS_UTEST_MPI_EPETRA_EXPONENTIAL_DATA_FIT_N5_HPP_

#include "Epetra_MpiComm.h"
#include "CONTAINERS_VECTOR"
#include "CONTAINERS_MATRIX"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.4
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f
 */

namespace pressio{ namespace solvers{ namespace test{

struct ExpDataFitN5
{
  const int numEq_ = 33;
  const int numUn_ = 5;

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
  shptr<nat_vec_type> tt_    = {}; // store the t value for the data
  shptr<nat_vec_type> yy_    = {}; // store the y values

  const std::vector<double> trueS = {0.37541005210628,
				     1.9358469126255,
				     -1.4646871365262,
				     0.012867534639885,
				     0.022122699662032};

  void storeTimes(){
    int i=0;
    for (auto const & it : myGel_){
      (*tt_)[i] = 10. * static_cast<scalar_type>(it);
      i++;
    };
  }

  void storeYValues(){
    const std::vector<double> allYs = {8.44e-1, 9.08e-1, 9.32e-1,
  				       9.36e-1, 9.25e-1, 9.08e-1,
  				       8.81e-1, 8.5e-1, 8.18e-1,
  				       7.84e-1, 7.51e-1, 7.18e-1,
  				       6.85e-1, 6.58e-1, 6.28e-1,
  				       6.03e-1, 5.8e-1, 5.58e-1,
  				       5.38e-1, 5.22e-1, 5.06e-1,
  				       4.9e-1, 4.78e-1, 4.67e-1,
  				       4.57e-1, 4.48e-1, 4.38e-1,
  				       4.31e-1, 4.24e-1, 4.2e-1,
  				       4.14e-1, 4.11e-1, 4.06e-1};

    int i=0;
    for (auto const & it : myGel_){
      (*yy_)[i] = allYs[it];
      i++;
    };
  }

  ExpDataFitN5(){
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
  }//setUp

  inline scalar_type model(const state_type & x, scalar_type t)const{
    return x[0] + x[1] * exp(-t*x[3]) + x[2]*exp(-t*x[4]);
  }

  void residual(const state_type& x, residual_type & R) const {
    for (auto i=0; i< NumMyElem_; i++){
      R[i] = (*yy_)[i] - this->model(x, (*tt_)[i]);
    };
  }

  residual_type residual(const state_type& x) const {
    residual(x, *R_);
    return *R_;
  }

  void jacobian(const state_type & x, jacobian_type & jac) const{
    for (int i=0; i<NumMyElem_; i++)
      {
	scalar_type t = (*tt_)[i];
	jac(i,0) = -1.0;
	jac(i,1) = -exp(-t*x[3]);
	jac(i,2) = -exp(-t*x[4]);
	jac(i,3) = x[1]*exp(-t*x[3])*t;
	jac(i,4) = x[2]*exp(-t*x[4])*t;
      }
    //jac.data()->Print(std::cout);
  }//end jacobian

  jacobian_type jacobian(const state_type & x) const{
    jacobian(x, *J_);
    return *J_;
  }//end jacobian

};

}}} //end namespace pressio::solvers::test
#endif




/* this is jacobian when using a crs matrix */
// void jacobian(const state_type & x, jacobian_type & jac) const{
//   std::vector<int> Indices {0, 1, 2, 3, 4};
//   std::vector<double> Values {0., 0., 0, 0, 0,};
//   for (int i=0; i<NumMyElem_; i++)
//     {
// 	scalar_type t = (*tt_)[i];

// 	int thisGID = myGel_[i];
// 	Values[0] = -1.0;
// 	Values[1] = -exp(-t*x[3]);
// 	Values[2] = -exp(-t*x[4]);
// 	Values[3] = x[1]*exp(-t*x[2])*t;
// 	Values[4] = x[2]*exp(-t*x[3])*t;

// 	if (jac.isFillingCompleted())
// 	  jac.data()->ReplaceGlobalValues(thisGID, 5, Values.data(), Indices.data() );
// 	else
// 	  jac.data()->InsertGlobalValues(thisGID, 5, Values.data(), Indices.data() );
//     }
//   if (!jac.isFillingCompleted())
//     jac.fillingIsCompleted();
// }//end jacobian
