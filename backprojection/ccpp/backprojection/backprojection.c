#include <backprojection.h>

#define KC (4 * M_PI / 3e8 * 5.8e9)
#define M_C 3e8
#define PLANS_FILENAME "fftw3Plans"
#define PROGRESS_STEP 10

//double phi2 =  (30 * M_PI / 180); // be careful here, this is the beamwidth divided by 2
#define PHI2 (1.482 / 2 * M_PI / 180) // be careful here, this is the beamwidth divided by 2

// backprojection using the analytical signal instead of the spectrum
// resampling calling resample2
int backProjection1(double* vec_x, int Nx,
                    double* vec_r, int Nr,
                    double* r_over, int Nover, double dx,
                    double complex * srf, int Naz, int Nf,
                    double* vec_az, double complex *img )
{
    double az;
    int ret=0;
    int naz;
    int xn;
    int rn;
    double d;
    int loop;
    int k;

    double complex x[Nf];
    double complex fftx[Nf];
    double complex y[Nover];
    double complex ffty[Nover];

    double complex aux1;
    double complex aux2;
    double complex aux4;

    fftw_plan px;
    fftw_plan py;

    px = fftw_plan_dft_1d(Nf, x, fftx, FFTW_FORWARD, FFTW_MEASURE);
    py = fftw_plan_dft_1d(Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

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
                if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(PHI2)) ) == 1. )
                {
                    d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                    aux1 = cexp( I * KC * d );
                    aux2 = myInterp( d, r_over, y, dx);
                    aux4 = aux1 * aux2;
                    img[xn * Nr + rn] += aux4;
                }
            }
        }
        loop++;
    }

    return ret;
}

// backprojection using the analytical signal instead of the spectrum
// resampling inline
int backProjection2(double* vec_x, int Nx,
                    double* vec_r, int Nr,
                    double* r_over, int Nover, double dx,
                    double complex* srf, int Naz, int Nf,
                    double* vec_az, double complex *img )
{
    double az;
    int ret=0;
    int naz;
    int xn;
    int rn;
    double d;
    int loop;

    double complex x[Nf];
    double complex fftx[Nf];
    double complex y[Nover];
    double complex ffty[Nover];

    double complex aux1;
    double complex aux2;
    double complex aux4;

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
                if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(PHI2)) ) == 1. )
                {
                    d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                    aux1 = cexp( I * KC * d );
                    aux2 = myInterp( d, r_over, y, dx);
                    aux4 = aux1 * aux2;
                    img[xn * Nr + rn] += aux4;
                }
            }
        }
        loop++;
    }

    return ret;
}

// backprojection using the analytical signal instead of the spectrum
// with OMP directives
// resampling calling resample2
int backProjectionOmp(double* vec_x, int Nx,
                      double* vec_r, int Nr,
                      double* r_over, int Nover, double dx,
                      double complex* srf, int Naz, int Nf,
                      double* vec_az, double complex *img )
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

    double complex *x = fftw_alloc_complex( numThreads * Nf );
    double complex *fftx = fftw_alloc_complex( numThreads * Nf );
    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan px[numThreads];
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
    {
        px[k] = fftw_plan_dft_1d(Nf, &x[k*Nf], &fftx[k*Nf], FFTW_FORWARD, FFTW_MEASURE);
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);
    }

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double az;
        int xn;
        int rn;
        double d;
        double complex aux1;
        double complex aux2;
        double complex aux4;
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
                    if ( pulse( (az-vec_x[xn]) / (vec_r[rn] * tan(PHI2)) ) == 1. )
                    {
                        d = sqrt( pow(vec_r[rn], 2.) + pow(az-vec_x[xn], 2.) );
                        aux1 = cexp( I * KC * d );
                        aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
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

int backProjectionOmpSlantRange(double* vec_az,
                                double* vec_rg,
                                double* r_over,
                                double complex* sr,
                                MyPosition *myPosition, double complex *img,
                                MyParameters params)
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

    int NAz = params.Nx;
    int NRg = params.Ny;
    int Nover = params.Nover;
    double dx = params.dx;
    int Naz = params.Naz;
    int Nf = params.Nf;
    double phi_a_deg = params.phi_a_deg;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double az;
        double rg;

        int azn;
        int rgn;
        double d = 0.;
        double daz;
        double valSin;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;

#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            az = myPosition[naz].x;
            rg = myPosition[naz].y;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            for (azn=0; azn<NAz; azn++)
            {
                daz = pow(az-vec_az[azn], 2.);
                for (rgn=0; rgn<NRg; rgn++)
                {
                    d = sqrt( daz + pow(rg-vec_rg[rgn], 2.) );
                    valSin = fabs( (az-vec_az[azn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[azn * NRg + rgn] += aux4;
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
                    printf("naz = %d (%.2f, %.2f) d = %.2f\n", naz, az, rg, d);
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

// back projection with ple zero padding
int backProjectionOmpGroundRange(double* vec_x,
                                 double* vec_r,
                                 double* r_over,
                                 double complex* sr,
                                 MyPosition *myPosition, double complex *img,
                                 MyParameters params)
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

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nover = params.Nover;
    double dx = params.dx;
    int Naz = params.Naz;
    int Nf = params.Nf;
    double hScene = params.hScene;
    double phi_a_deg = params.phi_a_deg;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d = 0.;
        double dxa;
        double dza;
        double valSin;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                for (rn=0; rn<Ny; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    valSin = fabs( (xa-vec_x[xn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Ny + rn] += aux4;
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

// this version takes the scene elevation as an argument
int backProjectionOmpGroundRange_LETG(double* vec_x,
                                      double* vec_y,
                                      double* vec_z,
                                      double* r_over,
                                      double complex* sr,
                                      MyPosition *myPosition, double complex *img,
                                      MyParameters_LETG params)
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;

    int Nx = params.Nx;
    int Ny = params.Ny;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double hScene = params.hScene;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;

    phi_max = phi_a_deg / 2. * M_PI / 180.;

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.2f\n", dx, kc );
    printf( "Nover = %d, Naz = %d\n", Nover, Naz );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n");
    printf( "schedule(static, 1)\n\n" );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int n;
        double d = 0.;
        //        double phi_az;
        double phi_s;
        double phi_a;
        double dxa;
        double dya;
        double dza;
        //        double valSin;
        double valCos;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(static, 1)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            for (n=0; n<(Nx*Ny); n++)
            {
                dxa = pow(xa-vec_x[n], 2.);
                dya = pow(ya-vec_y[n], 2.);
                dza = pow(za-vec_z[n], 2.);
                d = sqrt( dxa + dya + dza );
                //                valSin = fabs( (xa-vec_x[n]) / d );
                valCos = fabs( ( (vec_x[n] - xa) * params.uxx + (vec_y[n] - ya) * params.uxy ) ) / d;
                phi_s = acos(valCos);
                phi_a = M_PI/2 - phi_s;
                //                if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                if ( (phi_a < phi_max) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                {
                    aux1 = cexp( I * kc * d );
                    aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                    aux4 = aux1 * aux2;
#pragma omp critical
                    img[n] += aux4;
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    dxa = pow(xa-params.meanX, 2.);
                    dya = pow(ya-params.meanY, 2.);
                    dza = pow(za-hScene, 2.);
                    d = sqrt( dxa + dya + dza );
                    valCos = fabs( ( (params.meanX - xa) * params.uxx + (params.meanY - ya) * params.uxy ) ) / d;
                    phi_s = acos(valCos) * 180 / M_PI;
                    phi_a = 90 - phi_s;
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f, phi_s = %.2f, phi_a = %.2f 10h10\n", naz, xa, ya, za, d, phi_s, phi_a);
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

// this version takes the scene elevation as an argument and corrects for the motion during the pulse
int backProjectionOmpGroundRange_corr(double* scene_x,
                                      double* scene_y,
                                      double* scene_z,
                                      double* r_over,
                                      double complex* sr,
                                      MyPosition *myPosition, double* vel,
                                      double complex *img,
                                      MyParameters_LETG params)
{
    int ret=0;
    int naz;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;

    int Nx = params.Nx;
    int Ny = params.Ny;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double hScene = params.hScene;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;
    double lambdac;

    lambdac = 4 * M_PI / kc;
    phi_max = phi_a_deg / 2. * M_PI / 180.;

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.2f, lambdac = %.2f cm\n", dx, kc, lambdac * 100 );
    printf( "Nover = %d, Naz = %d\n", Nover, Naz );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n");
    printf( "schedule(static, 1)\n\n" );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        double vra;
        int n;
        double d = 0.;
        //        double phi_az;
        double phi_s;
        double phi_a;
        double dxa;
        double dya;
        double dza;
        double sinPhi;
        double cosPhi;
        double phi_corr;
        double t;
        double fd;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(static, 1)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;
            vra = vel[naz];

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            for (n=0; n<(Nx*Ny); n++)
            {
                dxa = pow(xa-scene_x[n], 2.);
                dya = pow(ya-scene_y[n], 2.);
                dza = pow(za-scene_z[n], 2.);
                d = sqrt( dxa + dya + dza );
                cosPhi = ( (scene_x[n] - xa) * params.uxx + (scene_y[n] - ya) * params.uxy ) / d;
                phi_s = acos(cosPhi);
                phi_a = M_PI/2-phi_s;
                sinPhi = sin(phi_a);
                fd = 2 * vra * sinPhi / lambdac; // Doppler shift
                t = 2 * d / M_C;
                phi_corr = - 2 * M_PI * fd * t;
                if ( (fabs(phi_a) < phi_max) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                {
                    aux1 = cexp( I * ( kc * d ) );
                    aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                    aux4 = aux1 * aux2;
#pragma omp critical
                    img[n] += aux4;
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    dxa = pow(xa-params.meanX, 2.);
                    dya = pow(ya-params.meanY, 2.);
                    dza = pow(za-hScene, 2.);
                    d = sqrt( dxa + dya + dza );
                    cosPhi = ( (params.meanX - xa) * params.uxx + (params.meanY - ya) * params.uxy ) / d;
                    phi_s = acos(cosPhi);
                    phi_a = M_PI/2-phi_s;
                    sinPhi = sin(phi_a);
                    fd = 2 * vra * sinPhi / lambdac; // Doppler shift
                    t = 2 * d / M_C;
                    phi_corr = - 2 * M_PI * fd * t;
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f, phi_s = %.2f, phi_a = %.2f\n",
                           naz, xa, ya, za, d, phi_s*180/M_PI, phi_a*180/M_PI);
                    printf("Vr = %.2f, Doppler shift = %.2f Hz, phi_corr = %.2f° NOT USED t0\n",
                           vra, fd, phi_corr*180/M_PI);
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

// this version takes the scene elevation as an argument
// the backprojection is performed on a volume
int backProjectionOmpGroundRange_PoSAR_GB(double* sceneX,
                                          double* sceneY,
                                          double* sceneZ,
                                          double* r_over,
                                          double complex* sr,
                                          MyPosition *myPosition, double complex *img,
                                          MyParametersPoSAR_GB params)
{
    int ret=0;
    unsigned int naz;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;
    double sin_phi2;

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nz = params.Nz;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );
    phi_max = phi_a_deg / 2. * M_PI / 180.;

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n*************************************\nbackProjectionOmpGroundRange_PoSAR_GB\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, kc );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double antX;
        double antY;
        double antZ;
        int n;
        double d = 0.;
        double phi_a;
        double dxa;
        double dya;
        double dza;
        double valSin;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            antX = myPosition[naz].x;
            antY = myPosition[naz].y;
            antZ = myPosition[naz].z;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            for (n=0; n<(Nx*Ny*Nz); n++)
            {
                dxa = pow(antX-sceneX[n], 2.);
                dya = pow(antY-sceneY[n], 2.);
                dza = pow(antZ-sceneZ[n], 2.);
                d = sqrt( dxa + dya + dza );
                valSin = fabs( (antX-sceneX[n]) / d );
                //                valCos = abs( ( (vec_x[n] - xa) * params.uxx + (vec_y[n] - ya) * params.uxy ) ) / d;
                //                phi_s = acos(valCos);
                //                phi_a = M_PI/2 - phi_s;
                if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    //                if ( (phi_a < phi_max) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                {
                    aux1 = cexp( I * kc * d );
                    aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                    aux4 = aux1 * aux2;
#pragma omp critical
                    img[n] += aux4;
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    dxa = pow(antX-sceneX[n], 2.);
                    dya = pow(antY-sceneY[n], 2.);
                    dza = pow(antZ-sceneZ[n], 2.);
                    d = sqrt( dxa + dya + dza );
                    valSin = fabs( (antX-sceneX[n]) / d );
                    phi_a = asin(valSin) * 180 / M_PI;
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f, phi_a = %.2f\n", naz, antX, antY, antZ, d, phi_a);
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

// this version takes the scene elevation as an argument
// the backprojection is performed on a volume, sceneX, Y and Z are vectors
int backProjectionOmpGroundRange_PoSAR_GBalt(double* sceneX,
                                             double* sceneY,
                                             double* sceneZ,
                                             double* r_over,
                                             double complex* sr,
                                             MyPosition *myPosition, double complex *img,
                                             MyParametersPoSAR_GB params)
{
    int ret=0;
    unsigned int naz;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;
    double sin_phi2;

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nz = params.Nz;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );
    phi_max = phi_a_deg / 2. * M_PI / 180.;

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n****************************************\nbackProjectionOmpGroundRange_PoSAR_GBalt\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, kc );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double antX;
        double antY;
        double antZ;
        int nx, ny, nz;
        double d = 0.;
        double phi_a;
        double dxa;
        double dya;
        double dza;
        double valSin;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            antX = myPosition[naz].x;
            antY = myPosition[naz].y;
            antZ = myPosition[naz].z;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            for (nx=0; nx<Nx; nx++){
                dxa = pow(antX-sceneX[nx], 2.);
                for (ny=0; ny<Ny; ny++){
                    dya = pow(antY-sceneY[ny], 2.);
                    for (nz=0; nz<Nz; nz++){
                        dza = pow(antZ-sceneZ[nz], 2.);
                        d = sqrt( dxa + dya + dza );
                        valSin = fabs( (antX-sceneX[nx]) / d );
                        //                valCos = abs( ( (vec_x[n] - xa) * params.uxx + (vec_y[n] - ya) * params.uxy ) ) / d;
                        //                phi_s = acos(valCos);
                        //                phi_a = M_PI/2 - phi_s;
                        if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                            //                if ( (phi_a < phi_max) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                        {
                            aux1 = cexp( I * kc * d );
                            aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                            aux4 = aux1 * aux2;
#pragma omp critical
                            img[nx * Ny * Nz + ny * Nz + nz] += aux4;
                        }
                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    dxa = pow(antX-sceneX[nx], 2.);
                    dya = pow(antY-sceneY[ny], 2.);
                    dza = pow(antZ-sceneZ[nz], 2.);
                    d = sqrt( dxa + dya + dza );
                    valSin = fabs( (antX-sceneX[nx]) / d );
                    phi_a = asin(valSin) * 180 / M_PI;
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) d = %.2f, phi_a = %.2f\n", naz, antX, antY, antZ, d, phi_a);
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

// this version takes the scene elevation as an argument
// the backprojection is performed on a volume, sceneX, Y and Z are vectors
// each thread has its own image
// each antenna Tx and Rx is processed using its real position
int backProjectionOmpGroundRange_PoSAR_GB_a(double *sceneX,
                                            double *sceneY,
                                            double *sceneZ,
                                            double *r_over,
                                            double complex *sr,
                                            MyPosition *positionRx,
                                            MyPosition *positionTx,
                                            double complex *img,
                                            MyParametersPoSAR_GB params)
{
    int ret=0;
    unsigned int naz=0;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;
    double sin_phi2;
    double meanx, meany, meanz;
    int kx, ky, kz;
    double d_sr, dRx_sr, dTx_sr;
    double yRx, zRx;
    double yTx, zTx;
    double dyRx, dzRx;
    double dyTx, dzTx;
    double dmin, dmax;

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nz = params.Nz;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );
    phi_max = phi_a_deg / 2. * M_PI / 180.;
    dmin = r_over[0] >= 0 ? r_over[0] : 0;
    dmax = r_over[Nover-1] ;

    maxThreads = omp_get_max_threads();
    //    numThreads = maxThreads / 2;
    numThreads = 8;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n****************************************\nbackProjectionOmpGroundRange_PoSAR_GBalt\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, kc );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    printf( "new Rev!\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

    // compute the center of the scene
    meanx = 0.;
    meany = 0.;
    meanz = 0.;
    for (k = 0; k < Nx; k++)
        meanx += sceneX[k];
    for (k = 0; k < Ny; k++)
        meany += sceneY[k];
    for (k = 0; k < Nz; k++)
        meanz += sceneZ[k];
    meanx /= Nx;
    meany /= Ny;
    meanz /= Nz;
    printf("mean point of the scene: ( %.2f, %.2f, %.2f )\n", meanx, meany, meanz);

#pragma omp parallel num_threads( numThreads )
    {
        double xRx, yRx, zRx;
        double xTx, yTx, zTx;
        double antAverX;
        int kx, ky, kz;
        double d;
        double dRx = 0.;
        double dTx = 0.;
        double phi_a;
        double dxRx, dyRx, dzRx;
        double dxTx, dyTx, dzTx;
        double valSin, valSinRx, valSinTx;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        double mAllocation;

        double complex *y;
        double complex *ffty;
        fftw_plan py;
        double complex *imgThread;
        int tid;

#pragma omp critical
        {
            // MEMORY ALLOCATIONS
            y    = fftw_alloc_complex( Nover );
            ffty = fftw_alloc_complex( Nover );
            py   = fftw_plan_dft_1d((int) Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);
            imgThread = (double complex *) calloc( Nx * Ny * Nz, sizeof(double complex) );

            tid = omp_get_thread_num();
            mAllocation = ( Nx * Ny * Nz * sizeof(double complex) ) / 1e6;
            printf("thread %d *** allocation of %.2f Mbytes\n", tid, mAllocation);
        }

#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xRx = positionRx[naz].x;
            yRx = positionRx[naz].y;
            zRx = positionRx[naz].z;
            xTx = positionTx[naz].x;
            yTx = positionTx[naz].y;
            zTx = positionTx[naz].z;
            antAverX = ( xRx + xTx ) / 2;

            zeroPaddingAndIfft_ple( py, &sr[naz*Nf], Nf, ffty, Nover);

            for ( kx = 0; kx < Nx; kx++ ){
                dxRx = pow( xRx - sceneX[kx], 2.);
                dxTx = pow( xTx - sceneX[kx], 2.);
                for ( ky = 0; ky < Ny; ky++ ){
                    dyRx = pow( yRx - sceneY[ky], 2.);
                    dyTx = pow( yTx - sceneY[ky], 2.);
                    for ( kz = 0; kz < Nz; kz++ ){
                        dzRx = pow( zRx - sceneZ[kz], 2.);
                        dzTx = pow( zTx - sceneZ[kz], 2.);
                        dRx = sqrt( dxRx + dyRx + dzRx );
                        dTx = sqrt( dxTx + dyTx + dzTx );
                        d = (dRx + dTx) / 2;
                        valSin   = fabs( ( antAverX - sceneX[kx] ) / d );
                        valSinRx = fabs( ( xRx - sceneX[kx] ) / dRx );
                        valSinTx = fabs( ( xTx - sceneX[kx] ) / dTx );
                        //                valCos = abs( ( (vec_x[n] - xa) * params.uxx + (vec_y[n] - ya) * params.uxy ) ) / d;
                        //                phi_s = acos(valCos);
                        //                phi_a = M_PI/2 - phi_s;
//                        if ( (valSin < sin_phi2) && (d >= dmin) && (d <= dmax) )
                            if ( (valSinRx < sin_phi2) && (valSinTx < sin_phi2) && (d >= dmin) && (d <= dmax) )
                            //                if ( (phi_a < phi_max) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                        {
                            aux1 = cexp( I * kc * d );
                            aux2 = myInterp( d, r_over, y, dx);
                            aux4 = aux1 * aux2;
                            imgThread[kx * Ny * Nz + ky * Nz + kz] += aux4;
                        }
                    }
                }
            }

            if (tid == 0)
            {
                steps_completed++;
                if (steps_completed > progressLevel)
                {
                    dxRx = pow( xRx - meanx, 2. );
                    dyRx = pow( yRx - meany, 2. );
                    dzRx = pow( zRx - meanz, 2. );
                    d = sqrt( dxRx + dyRx + dzRx );
                    valSin = fabs( (xRx-meanx) / d );
                    phi_a = asin(valSin) * 180 / M_PI;
                    printf( "%.0f%% \n", steps_completed / stepsPerThread * 100 );
                    progressLevel = progressLevel + progressStep;
                    printf("naz = %d (%.2f, %.2f, %.2f) dI = %.2f, phi_a = %.2f\n", naz, xRx, yRx, zRx, d, phi_a);
                }
            }
        }

#pragma omp critical
        {
            for (k=0; k<Nx*Ny*Nz;k++)
            {
                img[k] += imgThread[k];
            }
            // MEMORY DEALLOCATIONS
            printf("thread %d *** deallocation\n", tid);
            fftw_free( y );
            fftw_free( ffty );
            fftw_destroy_plan( py );
            free( imgThread );
        }
    }

    printf( "100%%\n" );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

// same as _a but only one thread as there is only on azimuth position
int backProjectionOmpGroundRange_PoSAR_GB_lha(double *sceneX,
                                              double *sceneY,
                                              double *sceneZ,
                                              double *r_over,
                                              double complex *sr,
                                              MyPosition *positionRx,
                                              MyPosition *positionTx,
                                              double complex *img,
                                              MyParametersPoSAR_GB params)
{
    int ret=0;
    unsigned int naz=0;
    int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double phi_max;
    double sin_phi2;
    double meanx, meany, meanz;
    double dmin, dmax;

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nz = params.Nz;
    unsigned int Nover = params.Nover;
    double dx = params.dx;
    unsigned int Naz = params.Naz;
    int Nf = params.Nf;
    double phi_a_deg = params.phi_a_deg;
    double kc = params.kc;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );
    phi_max = phi_a_deg / 2. * M_PI / 180.;
    dmin = r_over[0] >= 0 ? r_over[0] : 0;
    dmax = r_over[Nover-1] ;

    printf("*****************************************\nbackProjectionOmpGroundRange_PoSAR_GB_lha\n\n");
    printf( "dx = %f, kc = %.12f\n", dx, kc );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

    // compute the center of the scene
    meanx = 0.;
    meany = 0.;
    meanz = 0.;
    for (k = 0; k < Nx; k++)
        meanx += sceneX[k];
    for (k = 0; k < Ny; k++)
        meany += sceneY[k];
    for (k = 0; k < Nz; k++)
        meanz += sceneZ[k];
    meanx /= Nx;
    meany /= Ny;
    meanz /= Nz;

    printf("mean point of the scene: ( %.2f, %.2f, %.2f )\n", meanx, meany, meanz);

    double xRx, yRx, zRx;
    double xTx, yTx, zTx;
    double antAverX;
    int kx, ky, kz;
    double d;
    double dRx = 0.;
    double dTx = 0.;
    double dxRx, dyRx, dzRx;
    double dxTx, dyTx, dzTx;
    double valSin, valSinRx, valSinTx;
    double complex aux1;
    double complex aux2;
    double complex aux4;

    double complex *y;
    double complex *ffty;
    fftw_plan py;

    // MEMORY ALLOCATIONS
    y    = fftw_alloc_complex( Nover );
    ffty = fftw_alloc_complex( Nover );
    py   = fftw_plan_dft_1d((int) Nover, ffty, y, FFTW_BACKWARD, FFTW_MEASURE);

    naz = 0;
    xRx = positionRx[naz].x;
    yRx = positionRx[naz].y;
    zRx = positionRx[naz].z;
    xTx = positionTx[naz].x;
    yTx = positionTx[naz].y;
    zTx = positionTx[naz].z;
    printf("Rx (%.3f, %.3f, %.3f)\n", xRx, yRx, zRx);
    printf("Tx (%.3f, %.3f, %.3f)\n", xTx, yTx, zTx);
    antAverX = ( xRx + xTx ) / 2;

    zeroPaddingAndIfft_ple( py, &sr[naz*Nf], Nf, ffty, Nover);

    for ( kx = 0; kx < Nx; kx++ ){
        dxRx = pow( xRx - sceneX[kx], 2.);
        dxTx = pow( xTx - sceneX[kx], 2.);
        for ( ky = 0; ky < Ny; ky++ ){
            dyRx = pow( yRx - sceneY[ky], 2.);
            dyTx = pow( yTx - sceneY[ky], 2.);
            for ( kz = 0; kz < Nz; kz++ ){
                dzRx = pow( zRx - sceneZ[kz], 2.);
                dzTx = pow( zTx - sceneZ[kz], 2.);
                dRx = sqrt( dxRx + dyRx + dzRx );
                dTx = sqrt( dxTx + dyTx + dzTx );
                d = (dRx + dTx) / 2;
                valSin   = fabs( ( antAverX - sceneX[kx] ) / d );
                valSinRx = fabs( ( xRx - sceneX[kx] ) / dRx );
                valSinTx = fabs( ( xTx - sceneX[kx] ) / dTx );
//                if ( (valSinRx < sin_phi2) && (valSinTx < sin_phi2) && (d >= dmin) && (d <= dmax) )
//                {
                    aux1 = cexp( I * kc * d );
                    aux2 = myInterp( d, r_over, y, dx);
                    aux4 = aux1 * aux2;
                    img[kx * Ny * Nz + ky * Nz + kz] += aux4;
//                }
            }
        }
    }

    // MEMORY DEALLOCATIONS
    fftw_free( y );
    fftw_free( ffty );
    fftw_destroy_plan( py );

    printf( "100%%\n\n" );

    ret = fftw_export_wisdom_to_filename(PLANS_FILENAME);

    return ret;
}

int backProjectionOmpGroundRange_NED(double* vec_x,
                                     double* vec_r,
                                     double* r_over,
                                     double complex* sr,
                                     MyPosition *myPosition, double complex *img,
                                     MyParameters params,
                                     MyHeading *myHeading)
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

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nover = params.Nover;
    double dx = params.dx;
    int Naz = params.Naz;
    int Nf = params.Nf;
    double hScene = params.hScene;
    double phi_a_deg = params.phi_a_deg;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );

    maxThreads = omp_get_max_threads();
    numThreads = maxThreads / 2;
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );
    printf( "Very first run may be long due to the fftw plan calculation.\n\n");

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d(Nover, &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d = 0.;
        double dxa;
        double dza;
        double valSin;
        double scalarProduct;
        double scal1;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            zeroPaddingAndIfft_ple( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                scal1 = (xa-vec_x[xn]) * myHeading[naz].cos;
                for (rn=0; rn<Ny; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    //                    valSin = fabs( (xa-vec_x[xn]) / d );
                    // COMPUTE THE SCALAR PRODUCT, ASSUMPTION: THE PLANE IS HORIZONTAL
                    scalarProduct = scal1
                            + ((ya-vec_r[rn]) * myHeading[naz].sin);
                    valSin = fabs( sin( M_PI / 2. - acos( scalarProduct / d ) ) );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Ny + rn] += aux4;
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

// back projection with lff zero padding
int backProjectionOmpGroundRangeb(double* vec_x,
                                  double* vec_r,
                                  double* r_over,
                                  double complex* sr,
                                  MyPosition *myPosition, double complex *img,
                                  MyParameters params)
{
    int ret=0;
    int naz;
    unsigned int k;
    int maxThreads;
    unsigned int numThreads;
    double steps_completed = 0.;
    double progressLevel = 0.;
    double progressStep;
    double stepsPerThread;
    double sin_phi2;

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nover = params.Nover;
    double dx = params.dx;
    int Naz = params.Naz;
    int Nf = params.Nf;
    double hScene = params.hScene;
    double phi_a_deg = params.phi_a_deg;

    sin_phi2 = sin( phi_a_deg / 2. * M_PI / 180. );

    maxThreads = omp_get_max_threads();
    numThreads = (unsigned int) (maxThreads / 2);
    if (numThreads==0)
        numThreads = 1;
    printf("\n\n\n**************\nBACKPROJECTION\n\n");
    printf( "maxThreads = %d, numThreads = %d\n", maxThreads, numThreads );
    printf( "dx = %f, kc = %.12f\n", dx, KC );
    printf( "Nover = %d\n", Nover );
    printf( "d_min = %.2f, dmax = %.2f\n\n", r_over[0], r_over[Nover-1] );

    stepsPerThread = Naz / numThreads;
    progressStep = stepsPerThread / PROGRESS_STEP;

    double complex *y = fftw_alloc_complex( numThreads * Nover );
    double complex *ffty = fftw_alloc_complex( numThreads * Nover );
    fftw_plan py[numThreads];

    for (k=0; k<numThreads; k++)
        py[k] = fftw_plan_dft_1d((int) (Nover), &ffty[k*Nover], &y[k*Nover], FFTW_BACKWARD, FFTW_MEASURE);

    if (Nf%2!=0)
        printf("warning, Nf should be a multiple of 2\n");

#pragma omp parallel num_threads( numThreads )
    {
        double xa;
        double ya;
        double za;
        int xn;
        int rn;
        double d = 0.;
        double dxa;
        double dza;
        double valSin;
        double complex aux1;
        double complex aux2;
        double complex aux4;
        int tid;
#pragma omp for schedule(dynamic)
        for (naz=0; naz<Naz; naz++)
        {
            tid = omp_get_thread_num();
            xa = myPosition[naz].x;
            ya = myPosition[naz].y;
            za = myPosition[naz].z;

            zeroPaddingAndIfft_lff( py[tid], &sr[naz*Nf], Nf, &ffty[tid*Nover], Nover);

            dza = pow(za-hScene, 2.);

            for (xn=0; xn<Nx; xn++)
            {
                dxa = pow(xa-vec_x[xn], 2.);
                for (rn=0; rn<Ny; rn++)
                {
                    d = sqrt( dxa + pow(ya-vec_r[rn], 2.) + dza );
                    valSin = fabs( (xa-vec_x[xn]) / d );
                    if ( (valSin < sin_phi2) && (d >= r_over[0]) && (d <= r_over[Nover-1]) )
                    {
                        aux1 = cexp( I * KC * d );
                        aux2 = myInterp( d, r_over, &y[tid*Nover], dx);
                        aux4 = aux1 * aux2;
#pragma omp critical
                        img[xn * Ny + rn] += aux4;
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

// resampling creating fftw plans inside the function
// possibility to use the function in a Jupyter notebook
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

// resampling using fftw plans
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

// ple
int zeroPaddingAndIfft_ple( fftw_plan py,
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

// lff
int zeroPaddingAndIfft_lff( fftw_plan py,
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
        ffty[k+Ny-Nx] = fftx[k];
    }

    fftw_execute(py);

    return ret;
}

// ple
// zero padding and ifft creating fftw plans inside the function
int zeroPaddingAndIfft4( fftw_complex* fftx, int Nx,
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

// lff
// zero padding and ifft creating fftw plans inside the function
int zeroPaddingAndIfft4b( fftw_complex* fftx, int Nx,
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

double complex myInterp( double x, double *xp, double complex *fp, double dx )
{
    int idx1;
    int idx2;
    double complex y;

    idx1 = (int) (floor( (x-xp[0]) / dx ));
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

int importPlans(void)
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
