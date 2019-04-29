
#if defined HAVE_TRILINOS
#ifndef QR_TPETRA_BLOCK_MULTI_VECTOR_TSQR_IMPL_HPP_
#define QR_TPETRA_BLOCK_MULTI_VECTOR_TSQR_IMPL_HPP_

#include "../qr_rfactor_solve_impl.hpp"
#include "Tpetra_TsqrAdaptor.hpp"

namespace rompp{ namespace qr{ namespace impl{

template<
  typename matrix_t, typename R_t,
  int n, int m,
  typename MV_t,
  template<typename...> class Q_type>
class TpetraBlockMVTSQR<matrix_t, R_t, n, m, MV_t, Q_type, void>{

  using int_t	     = int;
  using sc_t	     = typename core::details::traits<matrix_t>::scalar_t;
  using lo_t	     = typename core::details::traits<matrix_t>::local_ordinal_t;
  using go_t	     = typename core::details::traits<matrix_t>::global_ordinal_t;
  using node_t	     = typename core::details::traits<matrix_t>::node_t;

  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;
  using Q_t	     = Q_type<MV_t>;

  using tpetra_mv_t = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  //using Q2_t	     = Q_type<tpetra_mv_t>;
  using tsqr_adaptor_type = Tpetra::TsqrAdaptor<tpetra_mv_t>;

public:
  TpetraBlockMVTSQR() = default;
  ~TpetraBlockMVTSQR() = default;

  void computeThinOutOfPlace(matrix_t & A) {
    auto nVecs	   = A.globalNumVectors();
    auto blockSize = A.data()->getBlockSize();
    createLocalRIfNeeded(nVecs);

    // this is the row map of the block MV
    auto & ArowMap = A.getDataMap();
    createQIfNeeded(ArowMap, blockSize, nVecs);

    // get the multivector
    auto mv = A.data()->getMultiVectorView();
    auto Qv = Qmat_->data()->getMultiVectorView();
    tsqrAdaptor_.factorExplicit(mv, Qv, *localR_.get(), false);
  }

  // void computeThinInPlace(matrix_t & A) {}

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
      qr::impl::solve<vector_t, trcp_mat, n>(rhs, this->localR_, y);
  }

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  		   vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  // if R_type != wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    ::rompp::mpl::enable_if_t<
  	      !core::meta::is_dense_matrix_wrapper_teuchos<T>::value and
	      !std::is_void<T>::value
  	      > * = nullptr>
  const T & getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  // if R_type == wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    ::rompp::mpl::enable_if_t<
  	      core::meta::is_dense_matrix_wrapper_teuchos<T>::value and
	      !std::is_void<T>::value
  	      > * = nullptr>
  const T & getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(*this->localR_, Teuchos::View);
    return *this->Rmat_;
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize){
    if (localR_.is_null() or
    	(localR_->numRows()!=newsize and localR_->numCols()!=newsize)){
      localR_ = Teuchos::rcp(new serden_mat_t(newsize, newsize) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int blockSize, int numVecs){
    if (!Qmat_ or !Qmat_->hasRowMapEqualTo(map) ){
      Qmat_ = std::make_shared<Q_t>(map, blockSize, numVecs);
    }
  }

private:
  tsqr_adaptor_type tsqrAdaptor_;
  trcp_mat localR_			= {};
  int computedRank_			= {};
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace rompp::qr::impl
#endif
#endif
