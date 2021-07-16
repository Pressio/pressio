
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_BlockVector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "pressio_containers.hpp"
#include <Tpetra_Map_decl.hpp>
#include "pressio_containers.hpp"

using vec_t = Tpetra::BlockVector<double>;
using w_t = pressio::containers::Vector<vec_t>;
using tcomm = Teuchos::Comm<int>;
using map_t = Tpetra::Map<>;


TEST(containers_vector_tpetra_block, Constructor1)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  vec_t a(*map, 3);
  a.putScalar(1.);

  w_t b(a);
  ASSERT_EQ( b.extent(0), 15 );
  // b should have same values of a
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );

  // check that b and a are different objects
  auto av = a.getVectorView().getLocalViewHost();
  auto bv = b.data()->getVectorView().getLocalViewHost();
  ASSERT_TRUE( av.use_count() != 0 );
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( bv.data() != nullptr );

  // changing b should not change a
  b.data()->putScalar(2.);
  ASSERT_EQ( b.data()->getVectorView().norm1(), 90. );
  ASSERT_EQ( a.getVectorView().norm1(), 45. );
}


TEST(containers_vector_tpetra_block, Constructor2)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  // create native tpetra vector
  vec_t a(*map, 3);
  a.putScalar(1.);
  auto aptr = a.getVectorView().getLocalViewHost().data();
  ASSERT_EQ( a.getVectorView().norm1(), 45. );

  // move into pressio wrapper
  w_t b(std::move(a));
  // b should have same values of a
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );
  // for block we don't have move, so ptr should not be same of a
  auto bptr = b.data()->getVectorView().getLocalViewHost().data();
  ASSERT_TRUE( bptr != aptr );

  // // explicitly call destructor JUST for testing purposes
  // // in this case this should work because of the move
  // a.~vec_t();

  // check that b is valid
  auto bv = b.data()->getVectorView().getLocalViewHost();
  ASSERT_TRUE( bv.use_count() != 0 );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra_block, Constructor3)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t b(*map, 3);
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
}


TEST(containers_vector_tpetra_block, CopyConstructor)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(*map,3);
  a.data()->putScalar(1.);
  ASSERT_EQ( a.extent(0), 15 );

  w_t b(a);
  ASSERT_EQ( b.extent(0), 15 );
  // b should have same values of a
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );
  // changing b should NOT change a
  b.data()->putScalar(2.);
  ASSERT_EQ( a.data()->getVectorView().norm1(), 45. );
  ASSERT_EQ( b.data()->getVectorView().norm1(), 90. );

  // for copy constr b and a are different things
  auto av = a.data()->getVectorView().getLocalViewHost();
  auto bv = b.data()->getVectorView().getLocalViewHost();
  ASSERT_TRUE( av.use_count() == bv.use_count() );
  ASSERT_TRUE( av.data() != bv.data() );
  ASSERT_TRUE( av.data() != nullptr );
  ASSERT_TRUE( bv.data() != nullptr );
}


TEST(containers_vector_tpetra_block, MoveConstructor)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t a(*map, 3);
  a.data()->putScalar(1.);
  auto aptr = a.data()->getVectorView().getLocalViewHost().data();

  w_t b(std::move(a));
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );
  // if we move, the ptr should be same of a but for tpetra block 
  // the move does nto work, so it uses copy
  auto bptr = b.data()->getVectorView().getLocalViewHost().data();
  ASSERT_TRUE( bptr != aptr );

  // check that values of b are still correct
  auto bv = b.data()->getVectorView().getLocalViewHost();
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
  ASSERT_TRUE( bv.data() != nullptr );
}

TEST(containers_vector_tpetra_block, MoveAssign)
{
  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<const map_t> map =
    Teuchos::rcp(new map_t(15, 0, comm));

  w_t b(*map, 3);
  b.data()->putScalar(1.);
  ASSERT_EQ( b.data()->getVectorView().norm1(), 45. );
  auto bptr1 = b.data()->getVectorView().getLocalViewHost().data();

  {
    w_t a(*map, 3);
    a.data()->putScalar(2.);
    ASSERT_EQ( a.data()->getVectorView().norm1(), 90. );
    auto aptr = a.data()->getVectorView().getLocalViewHost().data();

    // move assign
    b = std::move(a);
    ASSERT_EQ( b.data()->getVectorView().norm1(), 90. );
    // if we move, the ptr should be same of a but for block wrapper 
    // move is not implemented, it is a copy
    auto bptr2 = b.data()->getVectorView().getLocalViewHost().data();
    ASSERT_TRUE( bptr2 != aptr );
    ASSERT_FALSE( bptr1 != bptr2 );
  }

  auto bv = b.data()->getVectorView().getLocalViewHost();
  ASSERT_EQ( b.data()->getVectorView().norm1(), 90. );
  ASSERT_EQ( b.extent(0), 15 );
  ASSERT_EQ( b.extentLocal(0), 5 );
  ASSERT_TRUE( bv.data() != nullptr );
}


// TEST(containers_vector_tpetra_block, NativeMove)
// {
//   Teuchos::RCP<const tcomm> comm =
//     Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
//   Teuchos::RCP<const map_t> map =
//     Teuchos::rcp(new map_t(15, 0, comm));

//   vec_t a(*map, 3);
//   a.putScalar(1.);
//   ASSERT_EQ( a.getMap()->getGlobalNumElements(), 15 );
//   ASSERT_EQ( a.getVectorView().norm1(), 45. );

//   vec_t b(*map, 3);
//   b = a;
//   ASSERT_EQ( b.getMap()->getGlobalNumElements(), 15 );
//   ASSERT_EQ( b.getVectorView().norm1(), 45. );

//   // this does not work
//   vec_t c(std::move(a));
//   ASSERT_EQ( c.getMap()->getGlobalNumElements(), 15 );
//   ASSERT_EQ( c.getVectorView().norm1(), 45. );
// }
