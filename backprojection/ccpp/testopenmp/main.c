#include <stdio.h>
#include <omp.h>

#define N 1000
#define CHUNK_SIZE 100

int main(int argc, char *argv[])
{
    omp_set_dynamic(0);
    omp_set_num_threads(10);
    int i = 0, res = 0, tid;

#pragma omp parallel private(tid)
    {
#pragma omp for schedule(static, CHUNK_SIZE) reduction(+: res)
        for(i = 0; i<N; i++)
        {
            tid = omp_get_thread_num();
            res = res + i;
            printf("%d %d\n", tid, res);
        }
    }

    printf("res = %d\n", res);

    return 0;
}
