int sum(int a, int b) {
return a + b;
}

typedef struct {
   float real;
   float imag;
} complex ;

int backProjection( complex* in, complex* out, int size)
{
	int k;
	float xa;

	for (k=0; k<size; k++)
	{
		out[k].real = in[k].real;
		out[k].imag = in[k].imag;
	}

	for xa in xa_vec:
    if loop%1000 == 0:
        print( "{} / {}".format(loop, nbLoops ) )
    img += np.exp( 1j * kc * (r**2 + (xa-x)**2 )**0.5 ) \
    * np.interp( (r**2 + (xa-x)**2 )**0.5, r_over, signal.resample( srf[loop,:], nbPointsResampled  ) ) \
    * pulse( (xa-x) / (r*np.tan(phi)) )
    loop += 1

	return size;
}
