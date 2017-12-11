#ifndef BACKPROJECTION_H
#define BACKPROJECTION_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <omp.h>

typedef struct{
    double rampNumber;
    double timeStamp;
    double x;
    double y;
    double z;
} MyPosition;

typedef struct{
    int Nx;
    int Nr;
    int Nover;
    double dx;
    int Naz;
    int Nf;
    double hScene;
    double phi_a_deg;
} MyParameters;

int backProjection(double *vec_x, int Nx,
                   double *vec_r, int Nr,
                   double *r_over, int Nover, double dx,
                   complex *srf, int Naz, int Nf,
                   double *vec_az, complex* img);
int backProjectionOmp(double* vec_x, int Nx,
                      double* vec_r, int Nr,
                      double* r_over, int Nover, double dx,
                      complex* srf, int Naz, int Nf,
                      double* vec_az, complex *img );
int backProjectionOmp2(double* vec_x, int Nx,
                       double* vec_r, int Nr,
                       double* r_over, int Nover, double dx,
                       complex* srf, int Npos, int Nf,
                       MyPosition* myPosition, complex *img , double hScene);
int backProjection2(double *vec_x, int Nx,
                    double *vec_r, int Nr,
                    double *r_over, int Nover, double dx,
                    complex *srf, int Naz, int Nf,
                    double *vec_az, complex* img);

int backProjectionOmpGroundRange(double* vec_x, int Nx,
                                 double* vec_r, int Nr,
                                 double* r_over, int Nover, double dx,
                                 complex* sr, int Naz, int Nf,
                                 MyPosition *myPosition, complex *img,
                                 double hScene);
int call_backProjectionOmpGroundRange(double* vec_x,
                                      double* vec_r,
                                      double* r_over,
                                      complex* sr,
                                      MyPosition *myPosition, complex *img,
                                      MyParameters params);
int backProjectionOmpGroundRangeb(double* vec_x, int Nx,
                                  double* vec_r, int Nr,
                                  double* r_over, int Nover, double dx,
                                  complex* sr, int Naz, int Nf,
                                  MyPosition *myPosition, complex *img,
                                  double hScene);

int resample(fftw_complex *x, fftw_complex *fftx, int Nx, fftw_complex *y, fftw_complex *ffty, int Ny);
int resample2( fftw_plan px, fftw_plan py, fftw_complex* fftx, int Nx, fftw_complex* ffty, int Ny);
int resample3( fftw_plan py, fftw_complex* fftx, int Nx, fftw_complex* ffty, int Ny);
int resample3b( fftw_plan py, fftw_complex* fftx, int Nx, fftw_complex* ffty, int Ny);
int resample4(fftw_complex* fftx, int Nx, fftw_complex *y, fftw_complex* ffty, int Ny);

complex interp(double x, double *xp, complex *fp , double dx);
double pulse( double x );

int measureAndSavePlans(fftw_complex* x, fftw_complex* fftx, int Nx, fftw_complex* y, fftw_complex* ffty, int Ny);

int importPlans();
int fftwInitThreads();

#endif // BACKPROJECTION_H
