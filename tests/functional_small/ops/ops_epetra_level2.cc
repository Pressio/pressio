
#include "epetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

// convenient alias for nice test names
using ops_epetra = epetraMultiVectorGlobSize15Fixture;

template <typename T>
struct EpetraAdapter;

template <>
struct EpetraAdapter<std::shared_ptr<Epetra_Vector>> {
    EpetraAdapter(const std::shared_ptr<Epetra_Vector>& v): ptr_(v) {}
    auto& operator()(int i) { return (*ptr_)[i]; }

    std::shared_ptr<Epetra_Vector> ptr_;
};

template <>
struct EpetraAdapter<std::shared_ptr<Epetra_MultiVector>> {
    EpetraAdapter(const std::shared_ptr<Epetra_MultiVector>& v): ptr_(v) {}
    auto& operator()(int i, int j) { return (*ptr_)[j][i]; }

    std::shared_ptr<Epetra_MultiVector> ptr_;
};

template <typename T>
auto adapter(T &v) {
    return v;
}

template <>
auto adapter<std::shared_ptr<Epetra_Vector>>(std::shared_ptr<Epetra_Vector> &v) {
    return EpetraAdapter<std::shared_ptr<Epetra_Vector>>(v);
}

template <>
auto adapter<std::shared_ptr<Epetra_MultiVector>>(std::shared_ptr<Epetra_MultiVector> &v) {
    return EpetraAdapter<std::shared_ptr<Epetra_MultiVector>>(v);
}

template <typename T>
auto get_global_host_view(T &view, const Epetra_Import& /* importer */) {
  // deep copy is needed here, also for expressions
  // Note: Kokkos may be not initialized and Eigen may not be in build
  auto size = ::pressio::ops::extent(view, 0);
  Teuchos::SerialDenseVector<int, double> v_local(size);
  for (decltype(size) i = 0; i < size; ++i) {
    v_local(i) = view(i);
  }
  return v_local;
}

template <>
auto get_global_host_view<Epetra_MultiVector>(Epetra_MultiVector &mv, const Epetra_Import &importer) {
  auto mv_import = std::make_shared<Epetra_MultiVector>(importer.TargetMap(), mv.NumVectors());
  if (0 != mv_import->Import(mv, importer, Insert))
    throw std::runtime_error("Epetra import failed");
  return adapter(mv_import);
}

template <>
auto get_global_host_view<Epetra_Vector>(Epetra_Vector &v, const Epetra_Import &importer) {
  auto v_import = std::make_shared<Epetra_Vector>(importer.TargetMap(), v.NumVectors());
  if (0 != v_import->Import(v, importer, Insert))
    throw std::runtime_error("Epetra import failed");
  return adapter(v_import);
}

// adapter proxy ops
namespace pressio{ namespace ops{

template <typename T> size_t extent(const EpetraAdapter<T> &v, int i) { return ::pressio::ops::extent(*v.ptr_, i); }
template <typename T> void   scale(EpetraAdapter<T> &v, double a)     { ::pressio::ops::scale(*v.ptr_, a); }
template <typename T> void   set_zero(EpetraAdapter<T> &v)            { ::pressio::ops::set_zero(*v.ptr_); }

}}

// Note: we need definitions of adapter ops, becuase there are no delcarations like
//  <typename T> void set_zero(T &v);
// which could be used by shared routines.
#include "ops_shared_level2.hpp"

template <
    typename FixtureType,
    typename TransMode,
    typename AType,
    typename XType,
    typename YType,
    typename ScalarType
    >
void test_impl(FixtureType &test, TransMode trans, ScalarType alpha, AType A, XType x, ScalarType beta, YType &y) {
    // copy original values, fetch whole vector (from all ranks) to do the verification locally
    auto y_ref_h = get_global_host_view(y, *test.importer_);

    // call tested routine on device
    pressio::ops::product(trans, alpha, A, x, beta, y);

    // Note: fetch whole data (from all ranks) and run
    // the simplified reference routine locally
    auto A_h = get_global_host_view(A, *test.importer_);
    auto x_h = get_global_host_view(x, *test.importer_);
    auto y_h = get_global_host_view(y, *test.importer_);

    // Important: run vanilla_gemv() on all ranks, because it can execute
    // cross-rank operations like scale() or set_zero()
    vanilla_gemv(trans, alpha, A_h, x_h, beta, y_ref_h);

    if (test.rank_ == 0) {
        const size_t y_size = ::pressio::ops::extent(y_h, 0);
        for (size_t i = 0; i < y_size; ++i) {
          EXPECT_DOUBLE_EQ(y_h(i), y_ref_h(i));
        }
    }
}

template <
    typename FixtureType,
    typename TransMode,
    typename AType,
    typename XType,
    typename YType
    >
void test_impl(FixtureType &test, TransMode trans, AType A, XType x, YType y) {
    // non-zero alpha and beta
    test_impl(test, trans, 2., A, x, -1., y);
    // simulate beta=0 with uninitialized y containing NaN
    const auto nan = std::nan("0");
    ::pressio::ops::fill(y, nan);
    test_impl(test, trans, 2., A, x, 0., y);
    // simulate alpha=0 with NaN in input matrix
    auto a00 = A[0][0];
    A[0][0] = nan;
    test_impl(test, trans, 0., A, x, 0., y);
    // alpha=0, beta=1
    ::pressio::ops::fill(y, 2.);
    test_impl(test, trans, 0., A, x, 1., y);
    // restore original value
    A[0][0] = a00;
}

TEST_F(ops_epetra, mv_prod_teuchos_vector)
{
    Teuchos::SerialDenseVector<int, double> x_teuchos(numVecs_);
    x_teuchos = 1.;

    test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_teuchos, *y_epetra);
}

TEST_F(ops_epetra, mv_prod_storein_teuchos_vector)
{
    Teuchos::SerialDenseVector<int, double> y_teuchos(numVecs_);
    y_teuchos = 1.;

    test_impl(*this, ::pressio::transpose{}, *myMv_, *x_epetra, y_teuchos);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN

TEST_F(ops_epetra, mv_prod_eigen_vector)
{
    Eigen::VectorXd x_eigen(numVecs_);
    x_eigen.setConstant(1.);

    test_impl(*this, ::pressio::nontranspose{}, *myMv_, x_eigen, *y_epetra);
}

TEST_F(ops_epetra, mv_T_vector_storein_eigen_vector)
{
    Eigen::VectorXd y_eigen(numVecs_);
    y_eigen.setConstant(1.);

    test_impl(*this, ::pressio::transpose{}, *myMv_, *x_epetra, y_eigen);
}

#endif
