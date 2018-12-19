
#ifdef HAVE_TRILINOS
#ifndef QR_TPETRA_MULTI_VECTOR_HPP_
#define QR_TPETRA_MULTI_VECTOR_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"
#include "qr_rfactor_solve_impl.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include "AnasaziTsqrOrthoManager.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziTpetraAdapter.hpp"


namespace rompp{ namespace qr{

template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       R_type,
	       ::rompp::qr::TSQR,
	       Q_type,
	       typename
	       std::enable_if<
		 core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >::type
	       >
  : public QRSolverBase<QRSolver<matrix_type, R_type, ::rompp::qr::TSQR, Q_type>,
			Q_type<typename core::details::traits<matrix_type>::wrapped_t>,
			R_type,
			matrix_type>{

  using MV = typename core::details::traits<matrix_type>::wrapped_t;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using LO_t = typename core::details::traits<matrix_type>::local_ordinal_t;
  using GO_t = typename core::details::traits<matrix_type>::global_ordinal_t;
  using map_t = typename core::details::traits<matrix_type>::data_map_t;
  using node_t = typename core::details::traits<matrix_type>::node_t;
  using hexsp = typename core::details::traits<matrix_type>::host_exec_space_t;
  using Q_t = Q_type<MV>;

  using this_t = QRSolver<matrix_type, R_type, ::rompp::qr::TSQR, Q_type>;
  using base_t = QRSolverBase<this_t, Q_t, R_type,  matrix_type>;

  using base_t::Qmat_;
  using base_t::Rmat_;

  using MVTraits = Anasazi::MultiVecTraits<sc_t, MV>;
  using ortho_t = Anasazi::TsqrOrthoManager<sc_t, MV>;

  using mat_type = Teuchos::SerialDenseMatrix<int, sc_t>;
  using mat_ptr = Teuchos::RCP<mat_type>;

  const std::string label_ = "Anasazi";

public:
  QRSolver() {
    OM_ = std::make_shared<ortho_t>(label_);
  }

  ~QRSolver() = default;

private:

  template <typename vector_in_t,
	    typename vector_out_t,
	    core::meta::enable_if_t<
	      core::meta::is_core_vector_wrapper<vector_in_t>::value and
	      core::meta::is_core_vector_wrapper<vector_out_t>::value and
	      // the type vector in should be from same package as Q
	      core::details::traits<vector_in_t>::wrapped_package_identifier ==
		core::details::traits<Q_t>::wrapped_package_identifier
	      > * = nullptr
	    >
  void projectImpl(const vector_in_t & vecIn,
		   vector_out_t & vecOut) const{
    core::ops::dot( *Qmat_, vecIn, vecOut );
  }


  // non-const because Anasazi::TSQR can take a non-cost matrix
  // but here we do not modify it
  void computeThinImpl(matrix_type & A){
    // get number of cols
    auto n = A.globalNumVectors();
    // the row map
    auto ArowMap = A.getRCPDataMap();

    // create Ortho manager if not already existing
    if (!OM_)
      OM_ = std::make_shared<ortho_t>(label_);

    // create Q factor
    Qmat_ = std::make_shared<Q_t>(ArowMap, n);

    // create B
    localR_ = Teuchos::rcp(new mat_type(n,n) );

    const int initialX1Rank = OM_->normalizeOutOfPlace(*A.data(),
						       *Qmat_->data(),
						       localR_);
    Rmat_ = std::make_shared<R_type>(localR_->values());

    auto err = OM_->orthonormError(*Qmat_->data());
    std::cout << " Rank = "
	      << initialX1Rank << " "
	      <<  n << " "
	      << " error = " << err
	      << std::endl;
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y){
    ::rompp::qr::impl::solve(rhs, *Rmat_, y);
  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

private:
  friend base_t;
  std::shared_ptr< ortho_t > OM_ = nullptr;
  mat_ptr localR_                = {};

};//end class

}} // end namespace rompp::qr
#endif
#endif
