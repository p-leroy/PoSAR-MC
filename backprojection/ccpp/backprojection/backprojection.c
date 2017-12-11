#include <backprojection.h>

#define KC (4 * M_PI / 3e8 * 5.8e9)
#define PLANS_FILENAME "fftw3Plans"
#define PROGRESS_STEP 10

double phi2 =  (30 * M_PI / 180); // be careful here, this is the beamwidth divided by 2

int backProjection1(double* vec_x, int Nx,
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
    int k;

    complex x[Nf];
    complex fftx[Nf];
    complex y[Nover];
    complex ffty[Nover];

    complex aux1;
    complex aux2;
    complex aux4;

    fftw_plan px;
    fftw_plan py;

    px = fftw_plan_dft_1d(Nf, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    loop = 0;
    for (naz=0; naz<Naz; naz++)
    {

        if (loop%1000 == 0)
            printf( "%d / %d\n", loop, Naz );
        az = vec_az[naz];

        for(k=0; k<Nf; k++)
            x[k] = srf[naz*Nf+k];
        resample2( px, py, fftx, Nf, ffty, Nover);

        for (xn=0; xn<Nx; xn++)
        {
            for (rn=0; rn<Nr; rn++)
            {
                if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(phi2)) ) == 1. )
                {
                    d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                    aux1 = cexp( I * KC * d );
                    aux2 = interp( d, r_over, y, dx);
                    aux4 = aux1 * aux2;
                    img[xn * Nr + rn] += aux4;
                }
            }
        }
        loop++;
    }

    return ret;
}

int backProjectionOmp(double* vec_x, int Nx,
                      double* vec_r, int Nr,
                      double* r_over, int Nover, double dx,
                      complex* srf, int Naz, int Nf,
                      double* vec_az, complex *img )
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    complex *x = fftw_alloc_complex( numThreads * Nf );
    complex *fftx = fftw_alloc_complex( numThreads * Nf );
    complex *y = fftw_alloc_complex( numThreads * Nover );
    complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan px[numThreads];
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
    {
        px[k] = fftw_plan_dft_1d(Nf, &x[k*Nf], &fftx[k*Nf], FFTW_FORWARD, FFTW_MEASURE);
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);
    }

    if (Nf%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double az;
        int xn;
        int rn;
        double d;
        complex aux1;
        complex aux2;
        complex aux4;
        int tid;
#pragma omp for
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            az = vec_az[naz];

            for(k=0; k<Nf; k++)
                x[tid*Nf+k] = srf[naz*Nf+k];
            resample2( px[tid], py[tid], &fftx[tid*Nf], Nf, &ffty[tid*Nover], Nover);

            for (xn=0; xn<Nx; xn++)
            {
                for (rn=0; rn<Nr; rn++)
                {
                    if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(phi2)) ) == 1. )
                    {
                        d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                        aux1 = cexp( I * KC * d );
                        aux2 = interp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Nr + rn] += aux4;
                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                }
            }
        }
    }

    printf( "100%%\n" );

    fftw_free( x );
    fftw_free( fftx );
    fftw_free( y );
    fftw_free( ffty );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

// resample ple
int backProjectionOmpGroundRange(double* vec_x, int Nx,
                                 double* vec_r, int Nr,
                                 double* r_over, int Nover, double dx,
                                 complex* sr, int Naz, int Nf,
                                 MyPosition *myPosition, complex *img,
                                 double hScene)
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double sin_phi2;

    sin_phi2 = sin( phi2 );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    complex *y = fftw_alloc_complex( numThreads * Nover );
    complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d;
        double dxa;
        double dza;
        double valSin;
        complex aux1;
        complex aux2;
        complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        //        for (naz=55700; naz<55701; naz++)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            resample3( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            //            printf("vec_rd[750] = (%.12f, %.12f), y[40] = (%.12f, %.12f)\n",
            //                   creal(ffty[tid*Nover+750]), cimag(ffty[tid*Nover+750]),
            //                    creal(y[tid*Nover+40]), cimag(y[tid*Nover+40])
            //                    );

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                for (rn=0; rn<Nr; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    valSin = fabs( (xa-vec_x[xn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = interp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Nr + rn] += aux4;
                    }
                    //                    if ((rn==40) && (xn==25))
                    //                    {
                    //                        printf("d = %.12f, aux1 = (%.12f, %.12f), aux2 = (%.12f, %.12f)\n",
                    //                               d,
                    //                               creal(aux1), cimag(aux1),
                    //                               creal(aux2), cimag(aux2)
                    //                               );
                    //                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f\n", naz, xa, ya, za, d);
                }
            }
        }
    }

    printf( "100%%\n" );

    fftw_free( y );
    fftw_free( ffty );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

// call backProjectionOmpGroundRange with a set of parameters
int call_backProjectionOmpGroundRange(double* vec_x,
                                      double* vec_r,
                                      double* r_over,
                                      complex* sr,
                                      MyPosition *myPosition, complex *img,
                                      MyParameters params)
{
    phi2 = params.phi_a_deg * M_PI / 180. / 2.;
    printf( " 2 * phi2 = phi_a = %.2f\n", 2 * phi2 * 180 / M_PI );
    return backProjectionOmpGroundRange(vec_x, params.Nx,
                                 vec_r, params.Nr,
                                 r_over, params.Nover, params.dx,
                                 sr, params.Naz, params.Nf,
                                 myPosition, img,
                                 params.hScene);
}

// resample lff
int backProjectionOmpGroundRangeb(double* vec_x, int Nx,
                                  double* vec_r, int Nr,
                                  double* r_over, int Nover, double dx,
                                  complex* sr, int Naz, int Nf,
                                  MyPosition *myPosition, complex *img,
                                  double hScene)
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double sin_phi2;

    sin_phi2 = sin( phi2 );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    complex *y = fftw_alloc_complex( numThreads * Nover );
    complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d;
        double dxa;
        double dza;
        double valSin;
        complex aux1;
        complex aux2;
        complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        //        for (naz=55700; naz<55701; naz++)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            resample3b( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            //            printf("vec_rd[750] = (%.12f, %.12f), y[40] = (%.12f, %.12f)\n",
            //                   creal(ffty[tid*Nover+750]), cimag(ffty[tid*Nover+750]),
            //                    creal(y[tid*Nover+40]), cimag(y[tid*Nover+40])
            //                    );

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                for (rn=0; rn<Nr; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    valSin = fabs( (xa-vec_x[xn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = interp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Nr + rn] += aux4;
                    }
                    //                    if ((rn==40) && (xn==25))
                    //                    {
                    //                        printf("d = %.12f, aux1 = (%.12f, %.12f), aux2 = (%.12f, %.12f)\n",
                    //                               d,
                    //                               creal(aux1), cimag(aux1),
                    //                               creal(aux2), cimag(aux2)
                    //                               );
                    //                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f\n", naz, xa, ya, za, d);
                }
            }
        }
    }

    printf( "100%%\n" );

    fftw_free( y );
    fftw_free( ffty );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

// without interp
int backProjectionOmpGroundRange_v2(double* vec_x, int Nx,
                                    double* vec_r, int Nr,
                                    double* r_over, int Nover, double dx,
                                    complex* sr, int Naz, int Nf,
                                    MyPosition *myPosition, complex *img,
                                    double hScene)
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double sin_phi2;

    sin_phi2 = sin( phi2 );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    complex *y = fftw_alloc_complex( numThreads * Nover );
    complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d;
        double dxa;
        double dza;
        double valSin;
        double overDx;
        complex aux1;
        complex aux2;
        complex aux4;
        int tid;
        int idx;

        overDx = 1 / dx;

#pragma omp for schedule(dynamic)
        //        for (naz=55700; naz<55701; naz++)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            resample3( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            //            printf("vec_rd[750] = (%.12f, %.12f), y[40] = (%.12f, %.12f)\n",
            //                   creal(ffty[tid*Nover+750]), cimag(ffty[tid*Nover+750]),
            //                    creal(y[tid*Nover+40]), cimag(y[tid*Nover+40])
            //                    );

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                for (rn=0; rn<Nr; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    valSin = fabs( (xa-vec_x[xn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        idx = d * overDx;
                        aux2 = y[tid*Nover+idx];
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Nr + rn] += aux4;
                    }
                    //                    if ((rn==40) && (xn==25))
                    //                    {
                    //                        printf("d = %.12f, aux1 = (%.12f, %.12f), aux2 = (%.12f, %.12f)\n",
                    //                               d,
                    //                               creal(aux1), cimag(aux1),
                    //                               creal(aux2), cimag(aux2)
                    //                               );
                    //                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f\n", naz, xa, ya, za, d);
                }
            }
        }
    }

    printf( "100%%\n" );

    fftw_free( y );
    fftw_free( ffty );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

int resample(fftw_complex* x, fftw_complex* fftx, int Nx,
             fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int k;
    fftw_plan px;
    fftw_plan py;
    int ret;
    int fftxSplitIdx;;

    fftxSplitIdx = Nx / 2;

    ret = fftxSplitIdx;

    px = fftw_plan_dft_1d(Nx, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    fftw_execute(px);

    for (k=0; k<fftxSplitIdx; k++)
    {
        ffty[k] = fftx[k] / Nx;
    }
    for (k=fftxSplitIdx; k<Nx; k++)
    {
        ffty[Ny-Nx+k] = fftx[k] / Nx;
    }
    for (k=fftxSplitIdx; k<Ny-fftxSplitIdx; k++)
    {
        ffty[k] = 0;
    }

    fftw_execute(py);

    return ret;
}

int resample2( fftw_plan px, fftw_plan py,
               fftw_complex* fftx, int Nx,
               fftw_complex* ffty, int Ny)
{
    int k;
    int ret;

    ret = Ny;

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

int resample3( fftw_plan py,
               fftw_complex* fftx, int Nx,
               fftw_complex* ffty, int Ny)
{
    int k;
    int ret;

    ret = Ny;

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    for (k=0; k<Nx/2; k++)
    {
        ffty[k] = fftx[k];
    }
    for (k=Nx/2; k<Ny-Nx/2; k++)
    {
        ffty[k] = 0;
    }
    for (k=Nx/2; k<Nx; k++)
    {
        ffty[k+Ny-Nx] = fftx[k];

    }

    fftw_execute(py);

    return ret;
}

int resample3b( fftw_plan py,
                fftw_complex* fftx, int Nx,
                fftw_complex* ffty, int Ny)
{
    int k;
    int ret;
    int vec_ind;

    vec_ind = ceil( ( Nx + 1. ) / 2. );

    ret = Ny;

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    for (k=0; k<vec_ind; k++)
    {
        ffty[k] = fftx[k];
    }
    for (k=vec_ind; k<Ny-vec_ind; k++)
    {
        ffty[k] = 0;
    }
    for (k=vec_ind; k<Nx; k++)
    {
        ffty[Ny-Nx+k] = fftx[k];
    }

    fftw_execute(py);

    return ret;
}

int resample4( fftw_complex* fftx, int Nx,
               fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int k;
    int ret;

    fftw_plan py;

    ret = Nx/2;

    py = fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

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

int resample4b( fftw_complex* fftx, int Nx,
                fftw_complex* y, fftw_complex* ffty, int Ny)
{
    int k;
    int ret;
    int vec_ind;

    fftw_plan py;

    vec_ind = ceil( ( Nx + 1. ) / 2. );
    ret = vec_ind;

    py = fftw_plan_dft_1d(Ny, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    for (k=0; k<vec_ind; k++)
    {
        ffty[k] = fftx[k] / Nx;
    }
    for (k=vec_ind; k<Ny-vec_ind; k++)
    {
        ffty[k] = 0;
    }
    for (k=vec_ind; k<Nx; k++)
    {
        ffty[k+Ny-Nx] = fftx[k] / Nx;
    }

    fftw_execute(py);

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
    complex aux4;

    int k;
    fftw_plan px;
    fftw_plan py;

    px = fftw_plan_dft_1d(Nf, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nx%2!=0)
        printf("warning, Nx should be a multiple of 2\n");

    loop = 0;
    for (naz=0; naz<Naz; naz++)
    {

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

        for (xn=0; xn<Nx; xn++)
        {
            for (rn=0; rn<Nr; rn++)
            {
                if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(phi2)) ) == 1. )
                {
                    d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                    aux1 = cexp( I * KC * d );
                    aux2 = interp( d, r_over, y, dx);
                    aux4 = aux1 * aux2;
                    img[xn * Nr + rn] += aux4;
                }
            }
        }
        loop++;
    }

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
