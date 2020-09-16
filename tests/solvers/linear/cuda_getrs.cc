
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
//#include <cblas.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

TEST(solvers_linear_direct, getrs){

  constexpr int N = 8192;
  constexpr unsigned long BILLION = 1e9;
  struct timespec start ,stop;
  double accum;
  cublasStatus_t stat;
  cudaError cudaStatus;
  cusolverStatus_t cusolverStatus;
  cusolverDnHandle_t handle;
  double *A, *B1, *B;
  double *d_A, *d_B, *d_Work;
  int *d_pivot, *d_info, Lwork;
  int info_gpu = 0;
  A = (double*)malloc(N*N*sizeof(double));
  B = (double*)malloc(N*sizeof(double));
  B1 = (double*)malloc(N*sizeof(double));
  for(int i=0;i<N*N;i++) A[i]=rand()/(double)RAND_MAX;
  for(int i=0;i<N;i++) B[i] = 0.0;
  for(int i=0;i<N;i++) B1[i] = 1.0;
  double al=1.0, bet=0.0;
  int incx=1, incy=1;
  dgemv('N',N,N,al,A,N,B1,incx, bet,B,incy);
  cudaStatus = cudaGetDevice (0);
  cusolverStatus = cusolverDnCreate(&handle);
  cudaStatus = cudaMalloc((void**)&d_A,N*N*sizeof(double));
  cudaStatus = cudaMalloc((void**)&d_B, N*sizeof(double));

  cudaStatus = cudaMalloc((void**)&d_pivot, N*sizeof(int));
  cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
  cudaStatus = cudaMemcpy(d_A, A, N*N*sizeof(double), cudaMemcpyHostToDevice ); // copy
  cudaStatus = cudaMemcpy(d_B, B, N*sizeof(double), cudaMemcpyHostToDevice ); // copy
  cusolverStatus = cusolverDnDgetrf_bufferSize(handle , N, N, d_A, N, &Lwork); // compute buffer size and prep.memory
  cudaStatus=cudaMalloc((void**)&d_Work,Lwork*sizeof(double)); clock_gettime(CLOCK_REALTIME ,&start); // timer start

  // LU factorization of d_A, with partial pivoting and row // interchanges; row i is interchanged with row d_pivot(i);
  cusolverStatus = cusolverDnDgetrf(handle,N,N,d_A,N,d_Work, d_pivot, d_info);
  cusolverStatus = cusolverDnDgetrs(handle, CUBLAS_OP_N, N, 1, d_A, N, d_pivot, d_B,N, d_info);

  cudaStatus = cudaDeviceSynchronize ();
  clock_gettime(CLOCK_REALTIME ,&stop);

  accum=(stop.tv_sec-start.tv_sec)+(stop.tv_nsec-start.tv_nsec)/(double)BILLION;
  printf("getrf+getrs time: %lf sec.\n",accum);//pr.elaps.time

  cudaStatus = cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);

  printf("after getrf+getrs: info_gpu = %d\n", info_gpu);

  cudaStatus = cudaMemcpy(B1, d_B, N*sizeof(double), cudaMemcpyDeviceToHost);

  printf("solution: ");
  for (int i = 0; i < 5; i++) printf("%g, ", B1[i]);
  printf(" ..."); // print first components of the solution
  printf("\n");

  cudaStatus = cudaFree(d_A);
  cudaStatus = cudaFree(d_B);
  cudaStatus = cudaFree(d_pivot);
  cudaStatus = cudaFree(d_info);
  cudaStatus = cudaFree(d_Work);
  free(A); free(B); free(B1);
  cusolverStatus = cusolverDnDestroy(handle);
  cudaStatus = cudaDeviceReset ();
}
