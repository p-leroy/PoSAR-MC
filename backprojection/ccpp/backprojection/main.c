#include <stdio.h>
#include <stdlib.h>

#include<backprojection.h>

#define NB_POINTS 10
#define DX 0.1
#define DR 0.1
#define DAZ 0.1
#define N_X 10
#define N_R 5
#define N_AZ 10
#define N_F 1500
#define N_OVER 15000
#define H_SCENE 90

int main(void)
{
    printf("Hello World!\n");

    int complexSize;

    complexSize = sizeof(double) * 2;

    double *vec_x = malloc( N_X * sizeof(double) );
    double *vec_r = malloc( N_R * sizeof(double) );
    double *r_over = malloc( N_OVER * sizeof(double) );
    complex *sr = malloc( N_AZ * N_F * complexSize );
    MyPosition *positions = malloc( N_AZ * sizeof(MyPosition) );
    complex *img = malloc( N_X * N_R * complexSize );

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

    // initialize img
    for(int x=0; x<N_X; x++)
    {
        for(int r=0; r<N_R; r++)
        {
            img[x*N_R+r] = 0;
        }
    }

    printf("sizeof(complex) = %d\n", complexSize);

    backProjectionOmpGroundRange(vec_x, N_X,
                                 vec_r, N_R,
                                 r_over, N_OVER, DX,
                                 sr, N_AZ, N_F,
                                 positions, img,
                                 H_SCENE);

    free(vec_x);
    free(vec_r);
    free(r_over);
    free(sr);
    free(positions);
    free(img);

    return 0;
}
