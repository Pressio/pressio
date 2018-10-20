
#ifndef QR_HACKED_HPP_
#define QR_HACKED_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "../../CORE_ALL"
#ifdef HAVE_TRILINOS
#include <Epetra_Import.h>
#endif

namespace rompp{ namespace qr{ namespace hack{


#ifdef HAVE_TRILINOS
// overload for:
// the input data is a wrapper of an Epetra_MultiVector
// the Q factor is a Q_type wrapper of an Epetra_MultiVector
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem 
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type, Q_type, R_type,
	       typename
	       std::enable_if<
     core::meta::is_core_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_multi_vector_epetra<
		   typename core::details::traits<matrix_type>::wrapped_t
		   >::value and 
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem
		 >::type
	       >{
private:
  using MV = Epetra_MultiVector;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<MV>;

public:
  QRSolver() = default;
  ~QRSolver() = default;

  
  void compute(const matrix_type & A){
    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();

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

    //std::cout << *eA2W.data();
    // do QR in Eigen
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,m);
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

    // store Q into replicated Epetra_Multivector
    Q_t locQ(locMap,Qm.cols());
    for (int i=0;i<Qm.rows();i++)
      for (int j=0;j<Qm.cols();j++)
    	locQ(i,j) = Qm(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Qm.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->data()->Import(*locQ.data(), importer2, Insert);
    Qmat_->data()->Print(std::cout);

    // store R factor
    Rmat_ = std::make_shared<R_type>(Rm);
    
  }//end method
  //-----------------------------------------------------
  
  const Q_t & cRefQFactor() const {
    return *Qmat_;
  }//end method
  //-----------------------------------------------------
  
  const R_type & cRefRFactor() const {
    return *Rmat_;
  }//end method
  //-----------------------------------------------------
     
private:
  std::shared_ptr<Q_t> Qmat_;
  std::shared_ptr<R_type> Rmat_;

};//end class
#endif
      
      
}}} // end namespace rompp::qr::hack
#endif
