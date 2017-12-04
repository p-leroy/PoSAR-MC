#ifndef BACKPROJECTION_H
#define BACKPROJECTION_H

#include <complex.h>
#include <fftw3.h>
#include <stdio.h>

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
int backProjection2(double *vec_x, int Nx,
                    double *vec_r, int Nr,
                    double *r_over, int Nover, double dx,
                    complex *srf, int Naz, int Nf,
                    double *vec_az, complex* img);

int resample(fftw_complex *x, fftw_complex *fftx, int Nx, fftw_complex *y, fftw_complex *ffty, int Ny);
int resample2( fftw_plan px, fftw_plan py, fftw_complex* fftx, int Nx, fftw_complex* ffty, int Ny);
complex interp(double x, double *xp, complex *fp , double dx);
double pulse( double x );

int measureAndSavePlans(fftw_complex* x, fftw_complex* fftx, int Nx, fftw_complex* y, fftw_complex* ffty, int Ny);

int importPlans();
int fftwInitThreads();

#endif // BACKPROJECTION_H
