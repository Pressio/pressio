
#if defined HAVE_TRILINOS
#ifndef QR_EPETRA_MV_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_
#define QR_EPETRA_MV_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_

#include "../qr_rfactor_solve_impl.hpp"

namespace rompp{ namespace qr{ namespace impl{

/* partially specialize for when n and m are dynamic.
 * This just means that the R matrix is dynamic
*/
template<typename matrix_t, typename R_t,
	 typename MV_t, template<typename...> class Q_type>
class ModGramSchmidtMVEpetra<
  matrix_t, R_t, algebra::constants::dynamic,
  algebra::constants::dynamic, MV_t, Q_type, void>{

  using int_t	     = int;
  using sc_t	     = typename algebra::details::traits<matrix_t>::scalar_t;
  using eig_dyn_mat  =  Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using R_nat_t	     = algebra::Matrix<eig_dyn_mat>;
  using Q_t	     = Q_type<MV_t>;
  static constexpr sc_t one_ = static_cast<sc_t>(1);
  static constexpr sc_t zero_ = static_cast<sc_t>(0);

public:
  ModGramSchmidtMVEpetra() = default;
  ~ModGramSchmidtMVEpetra() = default;

  void computeThinOutOfPlace(matrix_t & A) {
    auto nVecs = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);

    sc_t rkkInv = zero_;
    for (auto k=0; k<A.globalNumVectors(); k++)
    {
      auto & ak = (*A.data())(k);
      ak->Norm2(&localR_(k,k));
      rkkInv = one_/localR_(k,k);

      auto & qk = (*Qmat_->data())(k);
      qk->Update( rkkInv, *ak, zero_ );

      for (auto j=k+1; j<A.globalNumVectors(); j++){
	auto & aj = (*A.data())(j);
	qk->Dot(*aj, &localR_(k,j));
	aj->Update(-localR_(k,j), *qk, one_);
      }
    }
  }

  // void computeThinInPlace(matrix_t & A) {
  //     // auto nVecs = A.globalNumVectors();
  //     // createLocalRIfNeeded(nVecs);
  //     // computedRank_ = OM_->normalize(*A.data(), localR_);
  //     // assert(computedRank_ == nVecs);
  // }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
    // auto vecSize = y.size();
    auto & Rm = localR_.data()->template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
	       vector_out_t & vecOut) const{
    algebra::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize){
    if (localR_.rows()!=newsize or localR_.cols()!=newsize){
      localR_ = R_nat_t(newsize, newsize);
      localR_.setZero();
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int cols){
    if (!Qmat_ or !Qmat_->hasRowMapEqualTo(map) )
      Qmat_ = std::make_shared<Q_t>(map, cols);
  }

private:
  R_nat_t localR_			= {};

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace rompp::qr::impl
#endif
#endif
