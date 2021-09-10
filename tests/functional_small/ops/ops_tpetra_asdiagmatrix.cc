
#include "tpetra_only_fixtures.hpp"
#include "pressio/ops.hpp"

// TEST(ops_tpetra, asDiagMatrixProductMultivector)
// {
//   //using NT = Tpetra::Vector<>::node_type;
//   using tcomm = Teuchos::Comm<int>;
//   using map_t = Tpetra::Map<>;
//   using vec_t = Tpetra::Vector<>;
//   using mv_t  = Tpetra::MultiVector<>;
//   using ST = typename vec_t::scalar_type;
//   // using LO = typename vec_t::local_ordinal_type;
//   // using GO = typename vec_t::global_ordinal_type;

//   Teuchos::RCP<const tcomm> comm = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
//   auto rank    = comm->getRank();
//   auto numProc = comm->getSize();
//   EXPECT_EQ(numProc,3);

//   int numLocalEntries = 2;
//   auto numGlobalEntries = numProc * numLocalEntries;
//   Teuchos::RCP<const map_t> map = Teuchos::rcp(new map_t(numGlobalEntries, 0, comm));

//   auto shift = rank*numLocalEntries;

//   // 1. create vector and fill with data
//   ::pressio::containers::Vector<vec_t> v(map);
//   auto vh = v.data()->getLocalViewHost();
//   std::vector<ST> v0(numGlobalEntries);
//   for (int i=0; i<6; ++i) v0[i]= (double) i;
//   for (int i=0; i<numLocalEntries; ++i) vh(i,0) = v0[shift+i];

//   // 2. create as_diagonal_matrix object from v
//   const auto vD = pressio::containers::as_diagonal_matrix(v);

//   // 3. create mv and fill with data
//   ::pressio::containers::MultiVector<mv_t> mv(map, 3);
//   auto mvh = mv.data()->getLocalViewHost();
//   Eigen::MatrixXd A(numGlobalEntries, 3);
//   A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
//   A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 1.;
//   A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.;
//   A(3,0) = 0.; A(3,1) = 1.; A(3,2) = 0.;
//   A(4,0) = 1.; A(4,1) = 0.; A(4,2) = 0.;
//   A(5,0) = 0.; A(5,1) = 1.; A(5,2) = 1.;

//   for (int i=0; i<numLocalEntries; ++i)
//     for (int j=0; j<3; ++j)
//       mvh(i,j) = A(shift+i,j);

//   // 4. compute product
//   ::pressio::containers::MultiVector<mv_t> C(map, 3);
//   pressio::ops::product(pressio::nontranspose(),
//       pressio::nontranspose(),
//       1., vD, mv, 0., C);

//   auto Ch = C.data()->getLocalViewHost();
//   for (int i=0; i<numLocalEntries; ++i)
//     for (int j=0; j<3; ++j)
//       std::cout << rank << " " << Ch(i,j) << std::endl;

//   // 5. check correctness
//   Eigen::MatrixXd G(numGlobalEntries, 3);
//   G(0,0) = 0.; G(0,1) = 0.; G(0,2) = 0.;
//   G(1,0) = 3.; G(1,1) = 2.; G(1,2) = 1.;
//   G(2,0) = 0.; G(2,1) = 0.; G(2,2) = 2.;
//   G(3,0) = 0.; G(3,1) = 3.; G(3,2) = 0.;
//   G(4,0) = 4.; G(4,1) = 0.; G(4,2) = 0.;
//   G(5,0) = 0.; G(5,1) = 5.; G(5,2) = 5.;

//   ASSERT_EQ( C.extent(0), 6 );
//   ASSERT_EQ( C.extent(1), 3 );
//   for (int i=0; i<numLocalEntries; ++i)
//     for (int j=0; j<3; ++j)
//       EXPECT_DOUBLE_EQ(Ch(i,j), G(shift+i, j));
// }