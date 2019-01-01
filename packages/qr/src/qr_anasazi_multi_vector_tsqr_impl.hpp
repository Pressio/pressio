
#if defined HAVE_TRILINOS and HAVE_ANASAZI_TSQR
#ifndef QR_ANASAZI_MULTI_VECTOR_TSQR_HPP_
#define QR_ANASAZI_MULTI_VECTOR_TSQR_HPP_

#include "qr_rfactor_solve_impl.hpp"

#include "AnasaziTsqrOrthoManager.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"


namespace rompp{ namespace qr{ namespace impl{

template<typename matrix_t, typename Q_t, typename R_t,
	 typename sc_t, typename MV>
class AnasaziMVTSQR{
protected:

  using this_t	     = AnasaziMVTSQR<matrix_t, Q_t, R_t, sc_t, MV>;
  using int_t	     = int;
  using ortho_t      = Anasazi::TsqrOrthoManager<sc_t, MV>;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;

  const std::string label_ = "Anasazi";

protected:
  AnasaziMVTSQR() = default;
  ~AnasaziMVTSQR() = default;

protected:

  void computeThinOutOfPlace(matrix_t & A) {
    auto nVecs = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);
    computedRank_ = OM_->normalizeOutOfPlace(*A.data(),
    					     *Qmat_->data(),
    					     localR_);
    assert(computedRank_ == nVecs);
  }

  void computeThinInPlace(matrix_t & A) {
    auto nVecs = A.globalNumVectors();
    createLocalRIfNeeded(nVecs);
    computedRank_ = OM_->normalize(*A.data(), localR_);
    assert(computedRank_ == nVecs);
  }

  template <typename vector_t, int n>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
    qr::impl::solve<vector_t, trcp_mat, n>(rhs, this->localR_, y);
  }


  // if R_type != wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    core::meta::enable_if_t<
  	      !core::meta::is_teuchos_serial_dense_matrix_wrapper<T>::value and
	      !std::is_void<T>::value
  	      > * = nullptr>
  const T & getCRefRFactor() const {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  // if R_type == wrapper of Teuchos::SerialDenseMatrix
  template <typename T = R_t,
  	    core::meta::enable_if_t<
  	      core::meta::is_teuchos_serial_dense_matrix_wrapper<T>::value and
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
  void createLocalRIfNeeded(int n){
    if (localR_.is_null() or
    	(localR_->numRows()!=n and localR_->numCols()!=n)){
      localR_ = Teuchos::rcp(new serden_mat_t(n, n) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int n){
    if (!Qmat_ or !Qmat_->hasRowMapEqualTo(map) )
      Qmat_ = std::make_shared<Q_t>(map, n);
  }

protected:
  std::shared_ptr< ortho_t > OM_	= std::make_shared<ortho_t>(label_);
  trcp_mat localR_			= {};
  int computedRank_			= {};

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace rompp::qr::impl
#endif
#endif
