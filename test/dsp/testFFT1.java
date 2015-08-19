package test.dsp;
import dsp.FFTProc;
import dsp.FFTProc3;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ijaux.Util;
import ijaux.datatype.Pair;
import ijaux.dsp.DSP;
import ijaux.dsp.FFT;
import static java.lang.Math.*;

public class testFFT1
 {
	// FFT by colls
	
	static float[] xc_re={18,    15,    11,    13, // nfft=4
		    -8,    -6,    -4,     7,
		     2,     5,     9,     9
		    -8,    -6,    -4,     7};
	
	static float[] xc_im={0,     0,     0,     0,
		    -8,    -5,    -1,    -2,
		     0,     0,     0,     0,
		     8,     5,     1,     2};
	/*static float[] xc_re={18,	15,	11,	13,  // nfft=3;
			-7.5f,	-4.5f,	-1,	7,
			-7.5f	-4.5f	-1,	7};
	
	
	static float[] xc_im={0,	0,	0,	0,
			 0.866025403784439f,	2.59807621135332f,	5.19615242270663f,	0,
			-0.866025403784439f,	-2.59807621135332f,	-5.19615242270663f,	0};*/
	
	// FFT by rows
	static float[] xr_re={15,    -2,    -7,    -2,
		    16,     7,     2,     7,
		    26,     2,     6,     2};
	
	static float[] xr_im={0,     7,     0,    -7,
		     0,    -3,     0,     3,
		     0,    -6,     0,     6};
	
	// 2D FFT
	
	static float[] x2_re={57,     7,     1,     7, // padding to 4x4
		    -11,    -7,   -13,    -1,
		     25,    -7,    -3,    -7,
		    -11,    -1,   -13,    -7};
	
	static float[] x2_im={0,    -2,     0,     2,
		    -16,     6,    -2,   -20,
		      0,     4,     0,    -4,
		     16,    20,     2,    -6};
	
/*	static float[] x2_re={57,    7,    1,    7, // no padding
		 	-6f,	-3.90192378864668f,	-11	-9.09807621135332f,
			-6f,	-9.09807621135332f,	-11	-3.90192378864668f};

	static float[] x2_im={0,	-2,	0,	2,
			 8.66025403784439f,	7.16987298107781f,	3.46410161513775f,	-15.8301270189222f,
			-8.66025403784439f,	15.8301270189222f,	-3.46410161513775f,	-7.16987298107781f};
*/

	
	static float[] arr={1,  2,  3,  9,
			8,  5,  1, 2,
			9 , 8 , 7,  2};
	
	static float[] x={1,	2,	3,	9,	8,	5,	1,	2};

	static float[] xr={31,	-14.0710678118655f,	5,	0.0710678118654755f,	
			-5,	0.0710678118654755f,	5,	-14.0710678118655f	};
	
	static float[] xi={0,	-4.82842712474619f,	4,	-0.828427124746190f,	
		0,	0.828427124746190f,	-4,	4.82842712474619f};
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new ImageJ();
		
		float[] row1={8,  5,  1, 2};
		float[] row2=row1.clone();
		int width=4;
		int height=3;
		FloatProcessor fp=new  FloatProcessor(width, height,arr);
		
		int nfft=FFT.nfft(row1.length);
		
		System.out.println ("\\ "+DSP.nextpow2(row1.length)+" nfft "+nfft);
		
 
		System.out.println ("********************************");
		System.out.println ("FFT");
		float[] imag= new float[row1.length];
		float[] imag2= new float[row1.length];
		Pair<float[], float[]> carr= FFTProc.fftC2C1d(row1, imag, -1, nfft);
		float[] re=carr.first;
		float[] im=carr.second;
		System.out.println ("\nReal part");
		Util.printFloatArray(re);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(im);
		carr= FFTProc.fftR2Cp1d(row2, -1, nfft);
		re=carr.first;
		 im=carr.second;
		System.out.println ("\nReal part");
		Util.printFloatArray(re);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(im);
		System.out.println ("********************************");
		System.out.println ("\nIFFT");
		nfft=FFT.nfft(xr.length);
		Pair<float[], float[]> carri= FFTProc.fftC2C1d(xr, xi, +1, nfft);
		
		System.out.println ("\nReal part");
		Util.printFloatArray(carri.first);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(carri.second);
		
	 	FFTProc3 proc =new FFTProc3(width*height);
		Pair<FloatProcessor, FloatProcessor> pfx=proc.fftR2C1D(fp, -1, 0); 
		width=pfx.first.getWidth();
		height=pfx.first.getHeight();
		
		
		ImageStack is1=new ImageStack(width, height);
		is1.addSlice("re", pfx.first);
		is1.addSlice("im", pfx.second);
		
		new ImagePlus("fft Ox", is1).show();
		
		Pair<FloatProcessor, FloatProcessor> pfy=proc.fftR2C1D(fp, -1, 1); 
		width=pfy.first.getWidth();
		height=pfy.first.getHeight();
		
		ImageStack is2=new ImageStack(width, height);
		is2.addSlice("re", pfy.first);
		is2.addSlice("im", pfy.second);
		
		new ImagePlus("fft Oy", is2).show();
		
		Pair<FloatProcessor, FloatProcessor> pf3=proc.fftR2C2D(fp, -1);
		//Pair<FloatProcessor, FloatProcessor> pf3=proc.fftC2C1D(pfy, -1, 0);
		width=pf3.first.getWidth();
		height=pf3.first.getHeight();
		ImageStack is=new ImageStack(width, height);
		is.addSlice("re", pf3.first);
		is.addSlice("im", pf3.second);
		
		new ImagePlus("fft 2D", is).show();

	}
	
	//-----------------------------------------------------------------------------
	// name: rfft()
	// desc: real value fft
	//
	//   these routines from the CARL software, spect.c
	//   check out the CARL CMusic distribution for more source code
	//
	//   if forward is true, rfft replaces 2*N real data points in x with N complex 
	//   values representing the positive frequency half of their Fourier spectrum,
	//   with x[1] replaced with the real part of the Nyquist frequency value.
	//
	//   if forward is false, rfft expects x to contain a positive frequency 
	//   spectrum arranged as before, and replaces it with 2*N real values.
	//
	//   N MUST be a power of 2.
	//
	//-----------------------------------------------------------------------------
	 boolean first = true ;
	void rfft( float[] x, int N, boolean forward )
	{
	   
	    float c1, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, temp, theta ;
	    float xr, xi ;
	    int i, i1, i2, i3, i4, N2p1 ;

	   
	    float   PI = (float) Math.PI;
	   
	     

	    theta = PI/N ;
	    wr = 1.f ;
	    wi = 0.f ;
	    c1 = 0.5f ;

	    if( forward )
	    {
	        c2 = -0.5f ;
	        cfft( x, N, forward ) ;
	        xr = x[0] ;
	        xi = x[1] ;
	    }
	    else
	    {
	        c2 = 0.5f ;
	        theta = -theta ;
	        xr = x[1] ;
	        xi = 0.f ;
	        x[1] = 0.f ;
	    }
	    double tmp=Math.sin( 0.5*theta );
	    wpr = (float) (-2.*tmp*tmp) ;
	    wpi = (float) Math.sin( theta ) ;
	    N2p1 = (N<<1) + 1 ;
	    
	    for( i = 0 ; i <= N>>1 ; i++ )
	    {
	        i1 = i<<1 ;
	        i2 = i1 + 1 ;
	        i3 = N2p1 - i2 ;
	        i4 = i3 + 1 ;
	        if( i == 0 )
	        {
	            h1r =  c1*(x[i1] + xr ) ;
	            h1i =  c1*(x[i2] - xi ) ;
	            h2r = -c2*(x[i2] + xi ) ;
	            h2i =  c2*(x[i1] - xr ) ;
	            x[i1] =  h1r + wr*h2r - wi*h2i ;
	            x[i2] =  h1i + wr*h2i + wi*h2r ;
	            xr =  h1r - wr*h2r + wi*h2i ;
	            xi = -h1i + wr*h2i + wi*h2r ;
	        }
	        else
	        {
	            h1r =  c1*(x[i1] + x[i3] ) ;
	            h1i =  c1*(x[i2] - x[i4] ) ;
	            h2r = -c2*(x[i2] + x[i4] ) ;
	            h2i =  c2*(x[i1] - x[i3] ) ;
	            x[i1] =  h1r + wr*h2r - wi*h2i ;
	            x[i2] =  h1i + wr*h2i + wi*h2r ;
	            x[i3] =  h1r - wr*h2r + wi*h2i ;
	            x[i4] = -h1i + wr*h2i + wi*h2r ;
	        }

	        wr = (temp = wr)*wpr - wi*wpi + wr ;
	        wi = wi*wpr + temp*wpi + wi ;
	    }

	    if( forward )
	        x[1] = xr ;
	    else
	        cfft( x, N, forward ) ;
	}




	//-----------------------------------------------------------------------------
	// name: cfft()
	// desc: complex value fft
	//
	//   these routines from CARL software, spect.c
	//   check out the CARL CMusic distribution for more software
	//
	//   cfft replaces float array x containing NC complex values (2*NC float 
	//   values alternating real, imagininary, etc.) by its Fourier transform 
	//   if forward is true, or by its inverse Fourier transform ifforward is 
	//   false, using a recursive Fast Fourier transform method due to 
	//   Danielson and Lanczos.
	//
	//   NC MUST be a power of 2.
	//
	//-----------------------------------------------------------------------------
	void cfft( float[] x, int NC, boolean forward )
	{
	    float wr, wi, wpr, wpi, theta, scale ;
	    int mmax, ND, m, i, j, delta ;
	    ND = NC<<1 ;
	    bit_reverse( x, ND ) ;
	    float   TWOPI = (float) (2.*Math.PI);
	    
	    for( mmax = 2 ; mmax < ND ; mmax = delta )
	    {
	        delta = mmax<<1 ;
	        theta = TWOPI/( forward? mmax : -mmax ) ;
	        double tmp= Math.sin( 0.5*theta );
	        wpr = (float) (-2.*tmp*tmp) ;
	        wpi = (float) Math.sin( theta ) ;
	        wr = 1.f ;
	        wi = 0.f ;

	        for( m = 0 ; m < mmax ; m += 2 )
	        {
	            float rtemp, itemp ;
	            for( i = m ; i < ND ; i += delta )
	            {
	                j = i + mmax ;
	                rtemp = wr*x[j] - wi*x[j+1] ;
	                itemp = wr*x[j+1] + wi*x[j] ;
	                x[j] = x[i] - rtemp ;
	                x[j+1] = x[i+1] - itemp ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;
	            }

	            wr = (rtemp = wr)*wpr - wi*wpi + wr ;
	            wi = wi*wpr + rtemp*wpi + wi ;
	        }
	    }

	    // scale output
	    scale = (float)(forward ? 1./ND : 2.) ;
	    for(int u=0; u< ND; u++ )
	    	x[u] *= scale ;
	    
	}




	//-----------------------------------------------------------------------------
	// name: bit_reverse()
	// desc: bitreverse places float array x containing N/2 complex values
//	       into bit-reversed order
	//-----------------------------------------------------------------------------
	void bit_reverse( float[] x, int N )
	{
	    float rtemp, itemp ;
	    int i, j, m ;
	    for( i = j = 0 ; i < N ; i += 2, j += m )
	    {
	        if( j > i )
	        {
	            rtemp = x[j] ; itemp = x[j+1] ; /* complex exchange */
	            x[j] = x[i] ; x[j+1] = x[i+1] ;
	            x[i] = rtemp ; x[i+1] = itemp ;
	        }

	        for( m = N>>1 ; m >= 2 && j >= m ; m >>= 1 )
	            j -= m ;
	    }
	}
	
	/* Function for generating Specialized sequence of twiddle factors */
	void tw_gen (double [] w, int n)
	{
		int i, j, k;
		double x_t, y_t, theta1, theta2, theta3;

		for (j = 1, k = 0; j <= n >> 2; j = j << 2)
		{
			for (i = 0; i < n >> 2; i += j)
			{
				theta1 = 2 * PI * i / n;
				x_t = cos (theta1);
				y_t = sin (theta1);
				w[k] = (float) x_t;
				w[k + 1] = (float) y_t;

				theta2 = 4 * PI * i / n;
				x_t = cos (theta2);
				y_t = sin (theta2);
				w[k + 2] = (float) x_t;
				w[k + 3] = (float) y_t;

				theta3 = 6 * PI * i / n;
				x_t = cos (theta3);
				y_t = sin (theta3);
				w[k + 4] = (float) x_t;
				w[k + 5] = (float) y_t;
				k += 6;
			}
		}
	}

	void split_gen (double[] pATable, double[] pBTable, int n){
		int i;

		for (i = 0; i < n; i++){
			pATable[2 * i] = 0.5 * (1.0 - sin (2 * PI / (double) (2 * n) * (double) i));
			pATable[2 * i + 1] = 0.5 * (-1.0 * cos (2 * PI / (double) (2 * n) * (double) i));
			pBTable[2 * i] = 0.5 * (1.0 + sin (2 * PI / (double) (2 * n) * (double) i));
			pBTable[2 * i + 1] = 0.5 * (1.0 * cos (2 * PI / (double) (2 * n) * (double) i));
		}
	}


  float[] FFT_Split (int n, float[]pIn, double[] pATable, double[] pBTable ){
	int i;
	float Tr, Ti;

	float[] pOut=new float[2*n];

	pIn[2 * n] = pIn[0];
	pIn[2 * n + 1] = pIn[1];


	for (i = 0; i < n; i++) {
		Tr = (float) (pIn[2 * i] * pATable[2 * i] 
				- pIn[2 * i + 1] * pATable[2 * i + 1] 
				+ pIn[2 * n - 2 * i] * pBTable[2 * i] 
				+ pIn[2 * n - 2 * i + 1] * pBTable[2 * i + 1]);
		Ti = (float) (pIn[2 * i + 1] * pATable[2 * i] 
				+ pIn[2 * i] * pATable[2 * i + 1] 
				+ pIn[2 * n - 2 * i] * pBTable[2 * i + 1] 
				- pIn[2 * n - 2 * i + 1] * pBTable[2 * i]);
		pOut[2 * i] = Tr;
		pOut[2 * i + 1] = Ti;
		// Use complex conjugate symmetry properties to get the rest..
		pOut[4 * n - 2 * i] = Tr;
		pOut[4 * n - 2 * i + 1] = -Ti;

	}
	pOut[2 * n] = pIn[0] - pIn[1];
	pOut[2 * n + 1] = 0;
	return pOut;
}

  float[] IFFT_Split (int n, float[] pIn, double[] pATable, double[] pBTable ){
	int i;
	float Tr, Ti;
	float[] pOut=new float[2*n];
	
	for (i = 0; i < n; i++)	{
		Tr = (float) (pIn[2 * i] * pATable[2 * i] + pIn[2 * i + 1] * pATable[2 * i + 1] + pIn[2 * n - 2 * i] * pBTable[2 * i] - pIn[2 * n - 2 * i + 1] * pBTable[2 * i + 1]);

		Ti = (float) (pIn[2 * i + 1] * pATable[2 * i] - pIn[2 * i] * pATable[2 * i + 1] - pIn[2 * n - 2 * i] * pBTable[2 * i + 1] - pIn[2 * n - 2 * i + 1] * pBTable[2 * i]);

		pOut[2 * i] = Tr;
		pOut[2 * i + 1] = Ti;
	}
	return pOut;
}
}
