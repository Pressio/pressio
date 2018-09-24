
#ifndef SVD_MULTI_VECTOR_EPETRA_HPP_
#define SVD_MULTI_VECTOR_EPETRA_HPP_

#include <memory>
#include "svd_solver_generic_base.hpp"
#include "../../CORE_ALL"
#include <Epetra_Import.h>

namespace rompp{ 
namespace svd{

template<typename matrix_type,
	 template <typename...> class lsv_type,
	 template <typename...> class rsv_type,
	 typename sval_type>
class Solver<matrix_type, lsv_type, rsv_type, sval_type,
	     typename
	     std::enable_if<
	       core::meta::is_multi_vector_epetra<
		 typename core::details::traits<matrix_type>::wrapped_t
		 >::value	       
	       >::type
	     >
  : public SolverBase< Solver<matrix_type, lsv_type,
			      rsv_type, sval_type> >
{

private:
  using MV = Epetra_MultiVector;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<MV>;
  using rsv_t = rsv_type<MV>;
  using sval_t = sval_type;

public:
  Solver() = default;
  ~Solver() = default;
  
private:

  template<svdType svd_enum_value,
	   typename std::enable_if<
	     svd_enum_value==svdType::truncated
	     >::type * = nullptr>
  void computeImpl(matrix_type & A, int t, sc_t tol){
    tol_ = tol;
    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();
    // truncated svd cannot have more than # of cols
    assert(t <= n);
    // number of left and right singular vectors
    nU_ = t;
    nV_ = t;

    // convert it to replicated eptra matrix 
    Epetra_LocalMap locMap(m, 0, A.commCRef());
    Epetra_Import importer(locMap, ArowMap);
    matrix_type A2(locMap, n);
    A2.data()->Import(*A.data(), importer, Insert);

    // store it into an Eigen matrix
    core::Matrix<Eigen::MatrixXd> eA2W(m,n);
    for (int i=0;i<m;i++)
      for (int j=0;j<n;j++)
    	eA2W(i,j) = A2(i,j);

    // do svd with eigen
    Eigen::BDCSVD<Eigen::MatrixXd> eSVD;
    eSVD.compute(*eA2W.data(), Eigen::ComputeFullU );//| Eigen::ComputeFullV );

    // store target t left singular vectors into replicated Epetra_Multivector
    lsv_t locU(locMap,t);
    for (int i=0;i<m;i++)
      for (int j=0;j<t;j++)
    	locU(i,j) = eSVD.matrixU()(i,j);

    // import from local to distributed
    lsv_ = std::make_shared<lsv_t>(ArowMap, t);
    Epetra_Import importer2(ArowMap, locMap);
    lsv_->data()->Import(*locU.data(), importer2, Insert);
    
  }//end method
  //-----------------------------------------------------
  
  const lsv_t & cRefLeftSingularVectorsImpl() const {
    return *lsv_;
  }//end method
  //-----------------------------------------------------
  
  const rsv_t & cRefRightSingularVectorsImpl() const {
    return *rsv_;
  }//end method
  //-----------------------------------------------------
  
  const sval_t & singularValuesImpl() const{
    return *sval_;
  }//end method
  //-----------------------------------------------------
   
private:
  friend SolverBase< Solver<matrix_type, lsv_type,
			    rsv_type, sval_type> >;

private:
  sc_t tol_;
  std::shared_ptr<lsv_t> lsv_;
  std::shared_ptr<rsv_t> rsv_;
  std::shared_ptr<sval_t> sval_;
  std::vector<sc_t> eigVals_;
  int nU_;
  int nV_;

};//end class
  
}//end namespace svd
}//end namespace rompp
#endif 
