#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(const int argc, const char* argv[])
{
  const int N=1024*1024*1024;
  const int M=4;
  double t1, t2;
  int checksum=0;

  printf("OpenMP threads: %d\n", omp_get_max_threads());

  //////////////////////////////////////////////////////////////////
  // Case 1: stack-allocated array
  t1=omp_get_wtime();
  checksum=0;
#pragma omp parallel
  { // Each openmp thread should have a private copy of
    // bins_thread_stack on the stack:
    int bins_thread_stack[M];
    for (int j=0; j<M; j++) bins_thread_stack[j]=0;
#pragma omp for
    for (int i=0; i<N; i++)
      { // Accumulating every M-th number in respective array element
        const int j=i%M;
        bins_thread_stack[j]++;
      }
#pragma omp critical
    for (int j=0; j<M; j++) checksum+=bins_thread_stack[j];
  }
  t2=omp_get_wtime();
  printf("Time with stack array: %12.3f sec, checksum=%d (must be %d).\n", t2-t1, checksum, N);
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Case 2: heap-allocated array
  t1=omp_get_wtime();
  checksum=0;
  #pragma omp parallel
  { // Each openmp thread should have a private copy of
    // bins_thread_heap on the heap:
    int* bins_thread_heap=(int*)malloc(sizeof(int)*M);
    for (int j=0; j<M; j++) bins_thread_heap[j]=0;
  #pragma omp for
    for (int i=0; i<N; i++)
      { // Accumulating every M-th number in respective array element
        const int j=i%M;
        bins_thread_heap[j]++;
      }
  #pragma omp critical
    for (int j=0; j<M; j++) checksum+=bins_thread_heap[j];
    free(bins_thread_heap);
  }
  t2=omp_get_wtime();
  printf("Time with heap  array: %12.3f sec, checksum=%d (must be %d).\n", t2-t1, checksum, N);
  //////////////////////////////////////////////////////////////////

  return 0;
}
