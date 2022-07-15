
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "pressio/rom_lspg.hpp"

TEST(rom_lspg, hypred_updater_vector)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  auto rank    = comm.MyPID();
  auto numProc = comm.NumProc();
  EXPECT_EQ(numProc,3);

  std::vector<double> vec_stencil_std =
    {44.,41,2.,5.,-1,-6.,45.,11.,4.,-10.,-8,2.2,3.5,6.,8.};

  std::vector<double> vec_sample_std = {1.,2.,3.,4.,5.,6.,7.};

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
  int my_gids_stencil[10];
  int my_size_stencil = 0;
  int shift_stencil = 0;

  int my_gids_sample[10];
  int my_size_sample = 0;
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

  Epetra_Map stencil_map(-1, my_size_stencil, my_gids_stencil, 0, comm);

  Epetra_Map sample_map(-1, my_size_sample, my_gids_sample, 0, comm);

  Epetra_Vector vec_stencil(stencil_map);
  for (int i=0; i<my_size_stencil; ++i){
    vec_stencil[i] = vec_stencil_std[shift_stencil+i];
  }

  Epetra_Vector vec_sample(sample_map);
  for (int i=0; i<my_size_sample; ++i){
    vec_sample[i] = vec_sample_std[shift_sample+i];
  }

  pressio::rom::lspg::impl::HypRedUpdaterTrilinos F;
  F.updateSampleMeshOperandWithStencilMeshOne(vec_sample,  2.,
						     vec_stencil, 3.);

  if (rank==0){
    EXPECT_DOUBLE_EQ(vec_sample[0], 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(vec_sample[1],  2.*3.+2.*2.);
  }

  else if (rank==1){
    EXPECT_DOUBLE_EQ(vec_sample[0], 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(vec_sample[1], 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(vec_sample[2], -8.*3.+5.*2);
  }

  else if (rank==2){
    EXPECT_DOUBLE_EQ(vec_sample[0], 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(vec_sample[1],  6.*3.+7.*2);
  }
}

TEST(rom_lspg, hypred_updater_multi_vector)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  auto rank    = comm.MyPID();
  auto numProc = comm.NumProc();
  EXPECT_EQ(numProc,3);

  std::vector<double> mv_stencil_std =
    {44.,41,2.,5.,-1,-6.,45.,11.,4.,-10.,-8,2.2,3.5,6.,8.};

  std::vector<double> mv_sample_std = {1.,2.,3.,4.,5.,6.,7.};

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
  int my_gids_stencil[10];
  int my_size_stencil = 0;
  int shift_stencil = 0;

  int my_gids_sample[10];
  int my_size_sample = 0;
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

  Epetra_Map stencil_map(-1, my_size_stencil, my_gids_stencil, 0, comm);
  Epetra_Map sample_map(-1, my_size_sample, my_gids_sample, 0, comm);

  Epetra_MultiVector mv_stencil(stencil_map, 4);
  for (int i=0; i<my_size_stencil; ++i){
    mv_stencil[0][i] = mv_stencil_std[shift_stencil+i];
    mv_stencil[1][i] = mv_stencil_std[shift_stencil+i]+1;
    mv_stencil[2][i] = mv_stencil_std[shift_stencil+i]+2;
    mv_stencil[3][i] = mv_stencil_std[shift_stencil+i]+3;
  }

  Epetra_MultiVector mv_sample(sample_map, 4);
  for (int i=0; i<my_size_sample; ++i){
    mv_sample[0][i] = mv_sample_std[shift_sample+i];
    mv_sample[1][i] = mv_sample_std[shift_sample+i]+1;
    mv_sample[2][i] = mv_sample_std[shift_sample+i]+2;
    mv_sample[3][i] = mv_sample_std[shift_sample+i]+3;
  }

  pressio::rom::lspg::impl::HypRedUpdaterTrilinos F;
  F.updateSampleMeshOperandWithStencilMeshOne(mv_sample,  2.,
						     mv_stencil, 3.);

  mv_stencil.Print(std::cout);
  if (rank==0){
    EXPECT_DOUBLE_EQ(mv_sample[0][0], 41.*3.+1.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[1][0], 42.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[2][0], 43.*3.+3.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[3][0], 44.*3.+4.*2.);

    EXPECT_DOUBLE_EQ(mv_sample[0][1],  2.*3.+2.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[1][1],  3.*3.+3.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[2][1],  4.*3.+4.*2.);
    EXPECT_DOUBLE_EQ(mv_sample[3][1],  5.*3.+5.*2.);
  }

  else if (rank==1){
    EXPECT_DOUBLE_EQ(mv_sample[0][0], 45.*3.+3.*2);
    EXPECT_DOUBLE_EQ(mv_sample[1][0], 46.*3.+4.*2);
    EXPECT_DOUBLE_EQ(mv_sample[2][0], 47.*3.+5.*2);
    EXPECT_DOUBLE_EQ(mv_sample[3][0], 48.*3.+6.*2);

    EXPECT_DOUBLE_EQ(mv_sample[0][1], 11.*3.+4.*2);
    EXPECT_DOUBLE_EQ(mv_sample[1][1], 12.*3.+5.*2);
    EXPECT_DOUBLE_EQ(mv_sample[2][1], 13.*3.+6.*2);
    EXPECT_DOUBLE_EQ(mv_sample[3][1], 14.*3.+7.*2);

    EXPECT_DOUBLE_EQ(mv_sample[0][2], -8.*3.+5.*2);
    EXPECT_DOUBLE_EQ(mv_sample[1][2], -7.*3.+6.*2);
    EXPECT_DOUBLE_EQ(mv_sample[2][2], -6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(mv_sample[3][2], -5.*3.+8.*2);
  }

  else if (rank==2){
    EXPECT_DOUBLE_EQ(mv_sample[0][0], 3.5*3.+6.*2);
    EXPECT_DOUBLE_EQ(mv_sample[1][0], 4.5*3.+7.*2);
    EXPECT_DOUBLE_EQ(mv_sample[2][0], 5.5*3.+8.*2);
    EXPECT_DOUBLE_EQ(mv_sample[3][0], 6.5*3.+9.*2);

    EXPECT_DOUBLE_EQ(mv_sample[0][1],  6.*3.+7.*2);
    EXPECT_DOUBLE_EQ(mv_sample[1][1],  7.*3.+8.*2);
    EXPECT_DOUBLE_EQ(mv_sample[2][1],  8.*3.+9.*2);
    EXPECT_DOUBLE_EQ(mv_sample[3][1],  9.*3.+10.*2);
  }

}
