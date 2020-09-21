
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "pressio_containers.hpp"
#include <Tpetra_Map_decl.hpp>
#include "pressio_containers.hpp"

TEST(containers_vector_tpetra, Constructor1)
{
  using vec_t = Tpetra::Vector<>;
  using w_t = pressio::containers::Vector<vec_t>;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  vec_t a(map);
  a.putScalar(1.);

  w_t b(a);
  ASSERT_EQ( b.data()->norm1(), 15. );

  auto av = a.getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() != 0 );
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, Constructor2)
{
  using vec_t = Tpetra::Vector<>;
  using w_t = pressio::containers::Vector<vec_t>;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  vec_t a(map);
  a.putScalar(1.);

  w_t b(std::move(a));
  ASSERT_EQ( b.data()->norm1(), 15. );

  auto av = a.getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() == 0 );
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( av.data() == bv.data() );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, CopyConstructor)
{
  using vec_t = Tpetra::Vector<>;
  using w_t = pressio::containers::Vector<vec_t>;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(map);
  a.data()->putScalar(1.);

  w_t b(a);
  ASSERT_EQ( b.data()->norm1(), 15. );

  auto av = a.data()->getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() == bv.use_count() );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( av.data() != nullptr );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, CopyAssign)
{
  using vec_t = Tpetra::Vector<>;
  using w_t = pressio::containers::Vector<vec_t>;

  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(map);
  a.data()->putScalar(1.);

  w_t b(map);
  b = a;
  ASSERT_EQ( b.data()->norm1(), 15. );

  auto av = a.data()->getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() == bv.use_count() );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( av.data() != nullptr );
  ASSERT_TRUE( bv.data() != nullptr );
}


// TEST(containers_vector_tpetra, MoveConstructor)
// {
//   using vec_t = Tpetra::Vector<>;
//   using w_t = pressio::containers::Vector<vec_t>;

//   using tcomm = Teuchos::Comm<int>;
//   using map_t = Tpetra::Map<>;

//   Teuchos::RCP<const tcomm> comm =
//     Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
//   Teuchos::RCP<const map_t> map =
//     Teuchos::rcp(new map_t(15, 0, comm));

//   w_t a(map);
//   auto av = a.data()->getLocalViewHost();
//   a.data()->putScalar(1.);
//   ASSERT_EQ( a.extent(0), 15 );
//   ASSERT_EQ( a.extentLocal(0), 5 );
//   ASSERT_TRUE( av.data() != nullptr );
// }


// TEST(containers_vector_tpetra, MoveAssign)
// {
//   using vec_t = Tpetra::Vector<>;
//   using w_t = pressio::containers::Vector<vec_t>;

//   using tcomm = Teuchos::Comm<int>;
//   using map_t = Tpetra::Map<>;

//   Teuchos::RCP<const tcomm> comm =
//     Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
//   Teuchos::RCP<const map_t> map =
//     Teuchos::rcp(new map_t(15, 0, comm));

//   w_t a(map);
//   a.data()->putScalar(1.);

//   w_t b(map);
//   b = std::move(a);
//   ASSERT_EQ( b.data()->norm1(), 15. );

//   // auto av = a.data()->getLocalViewHost();
//   // auto bv = b.data()->getLocalViewHost();
//   // ASSERT_TRUE( av.use_count() == bv.use_count() );
//   // //ASSERT_TRUE( av.data() != bv.data() );
//   // //ASSERT_TRUE( av.data() != nullptr );
//   // //ASSERT_TRUE( bv.data() != nullptr );
// }
