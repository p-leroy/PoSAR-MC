#include <backprojection.h>
#include <math.h>

#define KC (4 * M_PI / 3e8 * 5.8e9)
#define PLANS_FILENAME "fftw3Plans"
#define PHI (6 * M_PI / 180)

int backProjection(double* vec_x, int Nx,
                   double* vec_r, int Nr,
                   double* r_over, int Nover, double dx,
                   complex* srf, int Naz, int Nf,
                   double* vec_az, complex *img )
{
    double az;
    int ret=0;
    int naz;
    int x;
    int r;
    double d;
    int loop;

    complex fftx[Nf];
    complex y[Nover];
    complex ffty[Nover];

    complex aux1;
    complex aux2;
    double aux3;
    complex aux4;

    printf("\n\nbackProjection\n\n");

    printf("srf => (%.10f, %.10f) (%.10f, %.10f)\n",
           creal(srf[14995 * Nf]), cimag(srf[14995 * Nf]),
            creal(srf[14995 * Nf + 1]), cimag(srf[14995 * Nf + 1])
            );

    loop = 0;
    for (naz=14995; naz<15005; naz++)
    {
        //        printf("naz = %d, img[0, 0] = (%f, %f), ",naz,
        //               creal(img[0 * Nr + 0]), cimag(img[0 * Nr + 0]));
        printf("\nnaz = %d\n", naz);
        if (loop%1000 == 0)
            printf( "%d / %d\n", loop, Naz );
        az = vec_az[naz];
        resample( &srf[naz * Nf], fftx, Nf, y, ffty, Nover );
        for (x=0; x<Nx; x++)
        {
            for (r=0; r<Nr; r++)
            {
                d = sqrt( pow(vec_r[r], 2) + pow(az-vec_x[x], 2) );
                aux1 = cexp( I * KC * d );
                aux2 = interp( d, r_over, y, dx);
                aux3 = pulse( (az-vec_x[x]) / (vec_r[r] * tan(PHI)) );
                aux4 = aux1 * aux2 * aux3;
                img[x * Nr + r] += aux4;
                if ((r==20)&&(x==160))
                {
                    printf("d = %f, aux1 = (%.6f, %.6f), aux2 = (%.6f, %.6f)\n aux3 = %.6f, aux4 = (%.6f, %.6f)\n",
                           d,
                           creal(aux1), cimag(aux1),
                           creal(aux2), cimag(aux2),
                           aux3,
                           creal(aux4), cimag(aux4)
                           );
                    printf("srf => (%.10f, %.10f) (%.10f, %.10f)\n",
                           creal(srf[naz * Nf]), cimag(srf[naz * Nf]),
                            creal(srf[naz * Nf + 1]), cimag(srf[naz * Nf + 1])
                            );
                    printf("y => (%.10f, %.10f) (%.10f, %.10f)\n",
                           creal(y[1000]), cimag(y[1000]),
                            creal(y[1001]), cimag(y[1001])
                            );
                }
            }
        }
        loop++;
    }

    return ret;
}

int backProjection2(double* vec_x, int Nx,
                    double* vec_r, int Nr,
                    double* r_over, int Nover, double dx,
                    complex* srf, int Naz, int Nf,
                    double* vec_az, complex *img )
{
    double az;
    int ret=0;
    int naz;
    int xn;
    int rn;
    double d;
    int loop;

    complex x[Nf];
    complex fftx[Nf];
    complex y[Nover];
    complex ffty[Nover];

    complex aux1;
    complex aux2;
    double aux3;
    complex aux4;

    int k;
    fftw_plan px;
    fftw_plan py;

    px = fftw_plan_dft_1d(Nf, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    printf("\n\nbackProjection\n\n");

    //    printf("srf => (%.10f, %.10f) (%.10f, %.10f)\n",
    //           creal(srf[14995 * Nf]), cimag(srf[14995 * Nf]),
    //            creal(srf[14995 * Nf + 1]), cimag(srf[14995 * Nf + 1])
    //            );

    loop = 0;
    for (naz=0; naz<Naz; naz++)
    {
        //                printf("naz = %d, img[0, 0] = (%f, %f), ",naz,
        //                       creal(img[0 * Nr + 0]), cimag(img[0 * Nr + 0]));
        //        printf("\nnaz = %d\n", naz);
        if (loop%1000 == 0)
            printf( "%d / %d\n", loop, Naz );
        az = vec_az[naz];

        //***********
        //***********
        // RESAMPLING

        for (k=0; k<Nf; k++)
            x[k] = srf[naz * Nf + k];

        fftw_execute(px);

        for (k=0; k<Nf/2; k++)
        {
            ffty[k] = fftx[k] / Nf;
        }
        for (k=Nf/2; k<Nover-Nf/2; k++)
        {
            ffty[k] = 0;
        }
        for (k=Nf/2; k<Nf; k++)
        {
            ffty[k+Nover-Nf] = fftx[k] / Nf;
        }

        fftw_execute(py);

        //        printf("srf => (%.10f, %.10f) (%.10f, %.10f)\n",
        //               creal(srf[naz * Nf]), cimag(srf[naz * Nf]),
        //                creal(srf[naz * Nf + 1]), cimag(srf[naz * Nf + 1])
        //                );
        //        printf("x => (%.10f, %.10f) (%.10f, %.10f)\n",
        //               creal(x[0]), cimag(x[0]),
        //                creal(x[1]), cimag(x[1])
        //                );
        //        printf("y => (%.10f, %.10f) (%.10f, %.10f)\n",
        //               creal(y[1000]), cimag(y[1000]),
        //                creal(y[1001]), cimag(y[1001])
        //                );
        //
        //***********
        //***********

        for (xn=0; xn<Nx; xn++)
        {
            for (rn=0; rn<Nr; rn++)
            {
                d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                aux1 = cexp( I * KC * d );
                aux2 = interp( d, r_over, y, dx);
                aux3 = pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(PHI)) );
                aux4 = aux1 * aux2 * aux3;
                img[xn * Nr + rn] += aux4;
                //                img[xn * Nr + rn] += cexp( I * KC * d )
                //                        * interp( d, r_over, y, dx)
                //                        * pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(PHI)) );
                if( (xn == 160) && (rn == 20) )
                {
                    if (naz==Naz-1)
                    {
                        printf("d = %.16f, aux1 = (%.16f, %.16f), aux2 = (%.16f, %.16f)\n aux3 = %.16f, aux4 = (%.16f, %.16f)\n",
                               d,
                               creal(aux1), cimag(aux1),
                               creal(aux2), cimag(aux2),
                               aux3,
                               creal(aux4), cimag(aux4)
                               );
                    }
                }
            }
        }
        loop++;
    }

    return ret;
}

int resample(fftw_complex* x, fftw_complex* fftx, int Nx,
             fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int k;
    fftw_plan px;
    fftw_plan py;
    int ret;

    ret = Ny;

    px = fftw_plan_dft_1d(Nx, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    fftw_execute(px);

    for (k=0; k<Nx/2; k++)
    {
        ffty[k] = fftx[k] / Nx;
    }
    for (k=Nx/2; k<Ny-Nx/2; k++)
    {
        ffty[k] = 0;
    }
    for (k=Nx/2; k<Nx; k++)
    {
        ffty[k+Ny-Nx] = fftx[k] / Nx;
    }

    fftw_execute(py);

    return ret;
}

complex interp( double x, double *xp, complex *fp, double dx )
{
    int idx1;
    int idx2;
    complex y;

    idx1 = floor( (x-xp[0]) / dx );
    idx2 = idx1 + 1;

    y = (fp[idx2] - fp[idx1]) / dx * (x - xp[idx1]) + fp[idx1];

    return y;
}

double pulse( double x )
{
    int ret = 0.;

    if ( (-0.5<x) && (x<0.5) )
        ret = 1.;

    return ret;
}

int measureAndSavePlans(fftw_complex* x, fftw_complex* fftx, int Nx,
                        fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int ret;

    fftw_plan_dft_1d(Nx, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

int measurePlans(fftw_complex* x, fftw_complex* fftx, int Nx,
                 fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int ret=0;

    fftw_plan_dft_1d(Nx, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    return ret;
}

int importPlans()
{
    int ret;
    ret = fftw_import_wisdom_from_filename(PLANS_FILENAME);
    return ret;
}

int fftwInitThreads()
{
    //    return fftw_init_threads();
    return 0;
}
