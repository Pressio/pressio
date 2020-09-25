
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "pressio_containers.hpp"
#include <Tpetra_Map_decl.hpp>
#include "pressio_containers.hpp"

using vec_t = Tpetra::Vector<double>;
using w_t = pressio::containers::Vector<vec_t>;
using tcomm = Teuchos::Comm<int>;
using map_t = Tpetra::Map<>;


TEST(containers_vector_tpetra, Constructor1)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  vec_t a(map);
  a.putScalar(1.);

  w_t b(a);
  // b should have same values of a
  ASSERT_EQ( b.data()->norm1(), 15. );

  // check that b and a are different objects
  auto av = a.getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() != 0 );
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( bv.data() != nullptr );

  // changing b should not change a
  b.data()->putScalar(2.);
  ASSERT_EQ( b.data()->norm1(), 30. );
  ASSERT_EQ( a.norm1(), 15. );
}

TEST(containers_vector_tpetra, Constructor2)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  // create native tpetra vector
  vec_t a(map);
  a.putScalar(1.);
  auto aptr = a.getLocalViewHost().data();

  // move into pressio wrapper
  w_t b(std::move(a));
  // b should have same values of a
  ASSERT_EQ( b.data()->norm1(), 15. );
  // if we move, the ptr should be same of a
  auto bptr = b.data()->getLocalViewHost().data();
  ASSERT_TRUE( bptr == aptr );

  // explicitly call destructor JUST for testing purposes
  // in this case this should work because of the move
  a.~vec_t();

  // check that b is valid
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, Constructor3)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t b(map);
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
}

TEST(containers_vector_tpetra, CopyConstructor)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(map);
  a.data()->putScalar(1.);

  w_t b(a);
  // b should have same values of a
  ASSERT_EQ( b.data()->norm1(), 15. );
  // changing b should NOT change a
  b.data()->putScalar(2.);
  ASSERT_EQ( a.data()->norm1(), 15. );
  ASSERT_EQ( b.data()->norm1(), 30. );

  // for copy constr b and a are different things
  auto av = a.data()->getLocalViewHost();
  auto bv = b.data()->getLocalViewHost();
  ASSERT_TRUE( av.use_count() == bv.use_count() );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( av.data() != nullptr );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, MoveConstructor)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(map);
  a.data()->putScalar(1.);
  auto aptr = a.data()->getLocalViewHost().data();

  w_t b(std::move(a));
  ASSERT_EQ( b.data()->norm1(), 15. );
  // if we move, the ptr should be same of a
  auto bptr = b.data()->getLocalViewHost().data();
  ASSERT_TRUE( bptr == aptr );

  // explicitly call destructor JUST for testing purposes
  // in this case this should work because of the move
  a.~w_t();

  // check that values of b are still correct
  auto bv = b.data()->getLocalViewHost();
  ASSERT_EQ( b.data()->norm1(), 15. );
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra, MoveAssign)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t b(map);
  b.data()->putScalar(1.);
  ASSERT_EQ( b.data()->norm1(), 15. );
  auto bptr1 = b.data()->getLocalViewHost().data();

  {
    w_t a(map);
    a.data()->putScalar(2.);
    ASSERT_EQ( a.data()->norm1(), 30. );
    auto aptr = a.data()->getLocalViewHost().data();

    // move assign
    b = std::move(a);
    ASSERT_EQ( b.data()->norm1(), 30. );
    // if we move, the ptr should be same of a
    auto bptr2 = b.data()->getLocalViewHost().data();
    ASSERT_TRUE( bptr2 == aptr );
    ASSERT_FALSE( bptr1 == bptr2 );
  }

  auto bv = b.data()->getLocalViewHost();
  ASSERT_EQ( b.data()->norm1(), 30. );
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
  ASSERT_TRUE( bv.data() != nullptr );
}
