#include <stdio.h>

#include<backprojection.h>

#define NB_POINTS 10
#define DX 0.1
#define DR 0.1
#define DAZ 0.1
#define N_X 10
#define N_R 5
#define N_AZ 1

int main(int argc, char *argv[])
{
    printf("Hello World!\n");

    double vec_az[N_AZ];
    double vec_x[N_X];
    double vec_r[N_R];
    complex img[N_X][N_R];

    // initialize vec_x
    for(int k=0; k<N_X; k++)
    {
        vec_x[k] = k * DX;
    }

    // initialize vec_r
    for(int k=0; k<N_R; k++)
    {
        vec_r[k] = k * DR;
    }

    // initialize vec_az
    for(int k=0; k<N_AZ; k++)
    {
        vec_az[k] = k * DAZ;
    }

    // initialize img
    for(int x=0; x<N_X; x++)
    {
        for(int r=0; r<N_R; r++)
        {
            img[x][r] = 0;
        }
    }

    //backProjection(float *vec_x, int Nx, float *vec_r, int Nr, float *vec_az, int Naz, complex **out )
    backProjection( vec_x, N_X, vec_r, N_R, vec_az, N_AZ, img );

    printf( "%f, %f\n", creal(img[0][0]), cimag(img[0][0]) );

    return 0;
}
