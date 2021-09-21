
#include <gtest/gtest.h>
#include <Tpetra_BlockVector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>

#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, hypred_updater_vector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::BlockVector<>;
  using ST = typename vec_t::scalar_type;
  using LO = typename vec_t::local_ordinal_type;
  using GO = typename vec_t::global_ordinal_type;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank    = comm->getRank();
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  std::vector<ST> vec_stencil_std =
    {44.,41,2.,5.,-1,-6.,45.,11.,4.,-10.,-8,2.2,3.5,6.,8.};

  std::vector<ST> vec_sample_std = {1.,2.,3.,4.,5.,6.,7.};

  /*
     stencil mesh vector:
	     gids:
     rank=0: 0,1,2,3,4
     rank=1: 5,6,7,8,9,10
     rank=2: 11,12,13,14

     sample mesh vector:
	     gids:
     rank=0: 1,2
     rank=1: 6,7,10
     rank=2: 12,13
  */

  // create stencil mesh vector
  GO my_gids_stencil[10];
  LO my_size_stencil = 0;
  int shift_stencil = 0;

  GO my_gids_sample[10];
  LO my_size_sample = 0;
  int shift_sample = 0;

  if (rank==0){
    my_gids_stencil[0]=0;
    my_gids_stencil[1]=1;
    my_gids_stencil[2]=2;
    my_gids_stencil[3]=3;
    my_gids_stencil[4]=4;
    my_size_stencil = 5;
    shift_stencil = 0;

    my_gids_sample[0]=1;
    my_gids_sample[1]=2;
    shift_sample = 0;
    my_size_sample = 2;
  }
  else if(rank==1){
    my_gids_stencil[0]=5;
    my_gids_stencil[1]=6;
    my_gids_stencil[2]=7;
    my_gids_stencil[3]=8;
    my_gids_stencil[4]=9;
    my_gids_stencil[5]=10;
    my_size_stencil = 6;
    shift_stencil = 5;

    my_gids_sample[0]=6;
    my_gids_sample[1]=7;
    my_gids_sample[2]=10;
    shift_sample = 2;
    my_size_sample = 3;
  }
  else if(rank==2){
    my_gids_stencil[0]=11;
    my_gids_stencil[1]=12;
    my_gids_stencil[2]=13;
    my_gids_stencil[3]=14;
    my_size_stencil = 4;
    shift_stencil = 5+6;

    my_gids_sample[0]=12;
    my_gids_sample[1]=13;
    shift_sample = 5;
    my_size_sample = 2;
  }

  map_t stencil_map(-1, my_gids_stencil, my_size_stencil, 0, comm);
  map_t sample_map(-1, my_gids_sample, my_size_sample, 0, comm);

  int blocksize = 3;
  ST rowview[blocksize];

  vec_t vec_stencil(stencil_map, blocksize);
  for (LO i=0; i<my_size_stencil; ++i){
    rowview[0] = vec_stencil_std[shift_stencil+i];
    rowview[1] = rowview[0];
    rowview[2] = rowview[0];
    vec_stencil.replaceLocalValues(i, rowview);
  }

  vec_t vec_sample(sample_map, blocksize);
  for (LO i=0; i<my_size_sample; ++i){
    rowview[0] = vec_sample_std[shift_sample+i];
    rowview[1] = rowview[0];
    rowview[2] = rowview[0];
    vec_sample.replaceLocalValues(i, rowview);
  }

  pressio::rom::lspg::impl::HypRedUpdaterTrilinos F;
  F.updateSampleMeshOperandWithStencilMeshOne(vec_sample,  2.,
					      vec_stencil, 3.);

  if (rank==0){
    auto row0 = vec_sample.getLocalBlock(0);
    EXPECT_DOUBLE_EQ(row0(0), 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(row0(1), 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(row0(2), 41.*3.+1.*2.);

    auto row1 = vec_sample.getLocalBlock(1);
    EXPECT_DOUBLE_EQ(row1(0),  2.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(row1(1),  2.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(row1(2),  2.*3.+2.*2.);
  }

  else if (rank==1){
    auto row0 = vec_sample.getLocalBlock(0);
    EXPECT_DOUBLE_EQ(row0(0), 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(row0(1), 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(row0(2), 45.*3.+3.*2);

    auto row1 = vec_sample.getLocalBlock(1);
    EXPECT_DOUBLE_EQ(row1(0), 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(row1(1), 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(row1(2), 11.*3.+4.*2);

    auto row2 = vec_sample.getLocalBlock(2);
    EXPECT_DOUBLE_EQ(row2(0), -8.*3.+5.*2);
    EXPECT_DOUBLE_EQ(row2(1), -8.*3.+5.*2);
    EXPECT_DOUBLE_EQ(row2(2), -8.*3.+5.*2);
  }

  else if (rank==2){
    auto row0 = vec_sample.getLocalBlock(0);
    EXPECT_DOUBLE_EQ(row0(0), 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(row0(1), 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(row0(2), 3.5*3.+6.*2);

    auto row1 = vec_sample.getLocalBlock(1);
    EXPECT_DOUBLE_EQ(row1(0),  6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(row1(1),  6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(row1(2),  6.*3.+7.*2);
  }
}

TEST(rom_lspg, hypred_updater_multi_vector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using mv_t = Tpetra::BlockMultiVector<>;
  using ST = typename mv_t::scalar_type;
  using LO = typename mv_t::local_ordinal_type;
  using GO = typename mv_t::global_ordinal_type;

  Teuchos::RCP<const tcomm> comm =
    Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  auto rank    = comm->getRank();
  auto numProc = comm->getSize();
  EXPECT_EQ(numProc,3);

  std::vector<ST> mv_stencil_std =
    {44.,41,2.,5.,-1,-6.,45.,11.,4.,-10.,-8,2.2,3.5,6.,8.};

  std::vector<ST> mv_sample_std = {1.,2.,3.,4.,5.,6.,7.};

  /*
     stencil mesh mvtor:
	     gids:
     rank=0: 0,1,2,3,4
     rank=1: 5,6,7,8,9,10
     rank=2: 11,12,13,14

     sample mesh mvtor:
	     gids:
     rank=0: 1,2
     rank=1: 6,7,10
     rank=2: 12,13
  */

  // create stencil mesh mvtor
  GO my_gids_stencil[10];
  LO my_size_stencil = 0;
  int shift_stencil = 0;

  GO my_gids_sample[10];
  LO my_size_sample = 0;
  int shift_sample = 0;

  if (rank==0){
    my_gids_stencil[0]=0;
    my_gids_stencil[1]=1;
    my_gids_stencil[2]=2;
    my_gids_stencil[3]=3;
    my_gids_stencil[4]=4;
    my_size_stencil = 5;
    shift_stencil = 0;

    my_gids_sample[0]=1;
    my_gids_sample[1]=2;
    shift_sample = 0;
    my_size_sample = 2;
  }
  else if(rank==1){
    my_gids_stencil[0]=5;
    my_gids_stencil[1]=6;
    my_gids_stencil[2]=7;
    my_gids_stencil[3]=8;
    my_gids_stencil[4]=9;
    my_gids_stencil[5]=10;
    my_size_stencil = 6;
    shift_stencil = 5;

    my_gids_sample[0]=6;
    my_gids_sample[1]=7;
    my_gids_sample[2]=10;
    shift_sample = 2;
    my_size_sample = 3;
  }
  else if(rank==2){
    my_gids_stencil[0]=11;
    my_gids_stencil[1]=12;
    my_gids_stencil[2]=13;
    my_gids_stencil[3]=14;
    my_size_stencil = 4;
    shift_stencil = 5+6;

    my_gids_sample[0]=12;
    my_gids_sample[1]=13;
    shift_sample = 5;
    my_size_sample = 2;
  }

  (void)shift_stencil;
  (void)shift_sample;

  map_t stencil_map(-1, my_gids_stencil, my_size_stencil, 0, comm);
  map_t sample_map(-1, my_gids_sample, my_size_sample, 0, comm);
  int blocksize = 3;
  //ST rowview[blocksize];
  mv_t mv_stencil(stencil_map, blocksize, 4);
  mv_t mv_sample(sample_map, blocksize, 4);

  // finish it!

  pressio::rom::lspg::impl::HypRedUpdaterTrilinos F;
  F.updateSampleMeshOperandWithStencilMeshOne(mv_sample,  2.,
					      mv_stencil, 3.);
}
