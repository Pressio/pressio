
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>

#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, hypred_updater_vector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using vec_t = Tpetra::Vector<>;
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

  Teuchos::RCP<const map_t> stencil_map =
    Teuchos::rcp(new map_t(-1, my_gids_stencil, my_size_stencil, 0, comm));

  Teuchos::RCP<const map_t> sample_map =
    Teuchos::rcp(new map_t(-1, my_gids_sample, my_size_sample, 0, comm));

  const LO num_lel_stencil = stencil_map->getNodeNumElements();
  vec_t vec_stencil(stencil_map);
  auto vh = vec_stencil.getLocalViewHost();
  for (LO i=0; i<num_lel_stencil; ++i){
    vh(i,0) = vec_stencil_std[shift_stencil+i];
  }

  const LO num_lel_sample = sample_map->getNodeNumElements();
  vec_t vec_sample(sample_map);
  auto vh2 = vec_sample.getLocalViewHost();
  for (LO i=0; i<num_lel_sample; ++i){
    vh2(i,0) = vec_sample_std[shift_sample+i];
  }

  pressio::rom::lspg::impl::HypRedUpdater F;
  F.update_sample_mesh_operand_with_stencil_mesh_one(vec_sample,  2.,
						     vec_stencil, 3.);

  auto vh3 = vec_sample.getLocalViewHost();
  if (rank==0){
    EXPECT_DOUBLE_EQ(vh3(0,0), 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(vh3(1,0),  2.*3.+2.*2.);
  }

  else if (rank==1){
    EXPECT_DOUBLE_EQ(vh3(0,0), 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(vh3(1,0), 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(vh3(2,0), -8.*3.+5.*2);
  }

  else if (rank==2){
    EXPECT_DOUBLE_EQ(vh3(0,0), 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(vh3(1,0),  6.*3.+7.*2);
  }
}




TEST(rom_lspg, hypred_updater_multi_vector)
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using mv_t = Tpetra::MultiVector<>;
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

  Teuchos::RCP<const map_t> stencil_map =
    Teuchos::rcp(new map_t(-1, my_gids_stencil, my_size_stencil, 0, comm));

  Teuchos::RCP<const map_t> sample_map =
    Teuchos::rcp(new map_t(-1, my_gids_sample, my_size_sample, 0, comm));

  const LO num_lel_stencil = stencil_map->getNodeNumElements();
  mv_t mv_stencil(stencil_map, 4);
  auto vh = mv_stencil.getLocalViewHost();
  for (LO i=0; i<num_lel_stencil; ++i){
    vh(i,0) = mv_stencil_std[shift_stencil+i];
    vh(i,1) = mv_stencil_std[shift_stencil+i]+1;
    vh(i,2) = mv_stencil_std[shift_stencil+i]+2;
    vh(i,3) = mv_stencil_std[shift_stencil+i]+3;
  }

  const LO num_lel_sample = sample_map->getNodeNumElements();
  mv_t mv_sample(sample_map, 4);
  auto vh2 = mv_sample.getLocalViewHost();
  for (LO i=0; i<num_lel_sample; ++i){
    vh2(i,0) = mv_sample_std[shift_sample+i];
    vh2(i,1) = mv_sample_std[shift_sample+i]+1;
    vh2(i,2) = mv_sample_std[shift_sample+i]+2;
    vh2(i,3) = mv_sample_std[shift_sample+i]+3;
  }

  pressio::rom::lspg::impl::HypRedUpdater F;
  F.update_sample_mesh_operand_with_stencil_mesh_one(mv_sample,  2.,
						     mv_stencil, 3.);

  auto vh3 = mv_sample.getLocalViewHost();
  if (rank==0){
    EXPECT_DOUBLE_EQ(vh3(0,0), 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(vh3(0,1), 42.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(vh3(0,2), 43.*3.+3.*2.);
    EXPECT_DOUBLE_EQ(vh3(0,3), 44.*3.+4.*2.);

    EXPECT_DOUBLE_EQ(vh3(1,0),  2.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(vh3(1,1),  3.*3.+3.*2.);
    EXPECT_DOUBLE_EQ(vh3(1,2),  4.*3.+4.*2.);
    EXPECT_DOUBLE_EQ(vh3(1,3),  5.*3.+5.*2.);
  }

  else if (rank==1){
    EXPECT_DOUBLE_EQ(vh3(0,0), 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(vh3(0,1), 46.*3.+4.*2);
    EXPECT_DOUBLE_EQ(vh3(0,2), 47.*3.+5.*2);
    EXPECT_DOUBLE_EQ(vh3(0,3), 48.*3.+6.*2);

    EXPECT_DOUBLE_EQ(vh3(1,0), 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(vh3(1,1), 12.*3.+5.*2);
    EXPECT_DOUBLE_EQ(vh3(1,2), 13.*3.+6.*2);
    EXPECT_DOUBLE_EQ(vh3(1,3), 14.*3.+7.*2);

    EXPECT_DOUBLE_EQ(vh3(2,0), -8.*3.+5.*2);
    EXPECT_DOUBLE_EQ(vh3(2,1), -7.*3.+6.*2);
    EXPECT_DOUBLE_EQ(vh3(2,2), -6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(vh3(2,3), -5.*3.+8.*2);
  }

  else if (rank==2){
    EXPECT_DOUBLE_EQ(vh3(0,0), 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(vh3(0,1), 4.5*3.+7.*2);
    EXPECT_DOUBLE_EQ(vh3(0,2), 5.5*3.+8.*2);
    EXPECT_DOUBLE_EQ(vh3(0,3), 6.5*3.+9.*2);

    EXPECT_DOUBLE_EQ(vh3(1,0),  6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(vh3(1,1),  7.*3.+8.*2);
    EXPECT_DOUBLE_EQ(vh3(1,2),  8.*3.+9.*2);
    EXPECT_DOUBLE_EQ(vh3(1,3),  9.*3.+10.*2);
  }
}
