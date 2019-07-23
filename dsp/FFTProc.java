package dsp;


import ijaux.datatype.*;

import static dsp.DSP.*;
import static dsp.FFTUtil.*;
import static java.lang.Math.*;
 

/*
 *  FFT radix 2 implementations
 *  (C) Dimiter Prodanov, 
 */
public class FFTProc {
	
	private static final float tol=2E-15f;
	public final static double TWOPI=2.*Math.PI;
	
	double[] cosTable=null, sinTable=null, ctable=null;
	
	private int len=-1;
	
	public FFTProc(int n) {
		n=DSP.nfft(n);
		init(n, -1);
		len=n;		
	}

	/**
	 * @param n
	 */
	private void init(int n,  int sign) {
		Pair<double[], double[]> ptable=FFTUtil.expTable2(n, sign);
		cosTable=ptable.first;
		sinTable=ptable.second;
		ctable=FFTUtil.dPrecompute2(n);
	}
	
  
	public static boolean debug=false;
	
	/*
	private int[] nfftp(ImageProcessor ip) {
		int width=ip.getWidth();
		int height=ip.getHeight();
		
		//Util.printIntArray(dim);
		int nfftx=FFT.nfft(width);
		int nffty=FFT.nfft(height);
		
		if (debug)
		System.out.println (" nfftx "+nfftx+" nffty "+nffty);
		return new int[]{nfftx,nffty};
	}
	*/
	 
	/*
	 *  for compatibility with vanilla ImageJ
	 */
	
	public static Pair<float[], float[]> fftC2C1d(float[] real, float[] imag, final int sign, final int nfft) {	
		
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
			imag=zeroPaddEnd(imag);
		}
 		if (nfft<real.length) {
			real=FFTUtil.truncate(real,nfft);
			imag=FFTUtil.truncate(imag,nfft);
		} 
		if (real.length!= imag.length)
			throw new IllegalArgumentException("Array size mismatch " + real.length +" " +imag.length );

		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 
		
		// Bit reversal:
		final int len = real.length;
		//final int hlen = len/2;
		FFTUtil.bit_reverse(real);
		FFTUtil.bit_reverse(imag);
		
		double theta, tmp, alpha, beta, wr, wi;
		float tmpr, tmpi;
		//System.out.println("len "+ real.length);
		// Danielson-Lanczos algorithm:
		for (int N = 2, hN = 1; N <= len; hN = N, N <<= 1) {
			//System.out.println("N // "+ N+ " hN "+hN);
			theta = sign*TWOPI/N;
			//tmp = sin(0.5*theta);
			//alpha = -2.0*tmp*tmp;
			alpha =cos(theta)-1;
			beta = sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (int k=0; k<hN; k++) {
				//System.out.println("m // "+ k);				
				for (int i=k; i<len; i+=N) {
					int j=i+hN;
					//System.out.println("i "+ i +" j "+ j+" ; hN " +hN+" N "+ N);
					// butterfly					
					tmpr = (float) (wr*real[j] - wi*imag[j]);
					tmpi = (float) (wr*imag[j] + wi*real[j]);	
					real[j] =  (real[i] - tmpr);
					imag[j] =  (imag[i] - tmpi);
					real[i] += tmpr;
					imag[i] += tmpi;
				}
				tmp=wr;
				wr += tmp*alpha - wi*beta;
				wi += wi*alpha + tmp*beta;
			}
		}
		
		if (sign==1) {
			for (int i=0; i<len; i++) {
				real[i] /=(float)len;
				imag[i] /=(float)len;
				if (Math.abs(imag[i]) <=tol)
					imag[i]=0;
			}
		}
		
		return new Pair<float[], float[]>(real,imag);
	} // end
	
	/*
	 *  for compatibility with vanilla ImageJ
	 */
	
	public static Pair<float[], float[]> fftR2C1d(float[] real, final int sign, final int nfft) {
		
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
		}
 		if (nfft<real.length) {
			real=FFTUtil.truncate(real,nfft);
		} 
 		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 
 		float[] imag=new float[real.length];
 
		final int len = real.length;
		FFTUtil.bit_reverse(real);
		
		double tmp ;
		// Danielson-Lanczos algorithm:
		for (int N = 2, hN = 1; N <= len; hN = N, N <<= 1) {
			final double theta = sign*TWOPI/N;
			// tmp = sin(0.5*theta);
			//final double alpha = -2.0*tmp*tmp;
			final double alpha = cos(theta)-1;
			final double beta = sin(theta);
			double wr = 1.0;
			double wi = 0.0;
			for (int k=0; k<hN; k++) {
				for (int i=k, j=k+hN; i<len; i+=N, j+=N) {
					final float tmpr = (float) (wr*real[j] - wi*imag[j]);
					final float tmpi = (float) (wr*imag[j] + wi*real[j]);
					real[j] = (real[i] - tmpr);
					imag[j] = (imag[i] - tmpi);
					real[i] += tmpr;
					imag[i] += tmpi;				
				}
				tmp=wr;
				wr += tmp*alpha - wi*beta;
				wi += wi*alpha + tmp*beta;
			}
		}
		
		if (sign==1) {
			for (int i=0; i<len; i++) {
				real[i] /=(float)len;
				imag[i] /=(float)len;
				if (Math.abs(imag[i])<=tol)
					imag[i]=0;
			}
		}
		
		return new Pair<float[], float[]>(real,imag);
	} // end
		
 
	
	public static float[] twiddleAndUnfold (double [] ct, float[] X){
		int N=X.length;
		//System.out.println(N);
		//System.out.println("method 4");
		float[] G = new float[N<<1];
		for (int k = 0; k < N; k+=2){
			int u = N-k;
			if (1>k) u=u-N;	
			//System.out.println(k  +" "+u);
			//System.out.println(k+"  "+(k+1)+" "+ (u+1)+ " "+ (u)+" "+ (N-k-1));
			// Real part:
			G[k]= (float) (-ct[k+1]*X[u+1]-ct[k]*X[u]+X[u]-ct[k+1]*X[k+1]+ct[k]*X[k]);
			// Imaginary part:
			G[k+1]=(float)( ct[k  ]*X[u+1]-X[u+1]-ct[k+1]*X[u]+ct[k]*X[k+1]+X[k]*ct[k+1]);		
			//System.out.println("// "+G[k  ]+"  "+G[k+1]);
		}		
		//last real
		G[N]=X[0]-X[1];
		//last imaginary
		G[N+1]=0;
		// unfolding the full FFT
		for (int i = 0; i<N; i+=2){
			int  k = N*2-1-i;// full measure
			//System.out.println(i+" "+(k));
			// Real part:
			G[k-1] = G[2+i];
			// Imaginary part:
			G[k] = - G[2+i+1]; 
		}
		return G; 
	}
	
	
	public static double[] twiddleAndUnfold(double [] ct, double[] g2){
		int N=g2.length;
		//System.out.println(N);
		//System.out.println("method 4");
		double[] G = new double[N<<1];
		for (int k = 0; k < N; k+=2){
			int u = N-k;
			if (1>k) u=u-N;	
			//System.out.println(k  +" "+u);
			//System.out.println(k+"  "+(k+1)+" "+ (u+1)+ " "+ (u)+" "+ (N-k-1));
			// Real part:
			G[k]=(-ct[k+1]*g2[u+1]-ct[k]*g2[u]+g2[u]-ct[k+1]*g2[k+1]+ct[k]*g2[k]);
			// Imaginary part:
			G[k+1]=( ct[k  ]*g2[u+1]-g2[u+1]-ct[k+1]*g2[u]+ct[k]*g2[k+1]+g2[k]*ct[k+1]);		
			//System.out.println("// "+G[k  ]+"  "+G[k+1]);
		}		
		//last real
		G[N]=g2[0]-g2[1];
		//last imaginary
		G[N+1]=0;
		// unfolding the full FFT
		for (int i = 0; i<N; i+=2){
			int  k = N*2-1-i;// full measure
			//System.out.println(i+" "+(k));
			// Real part:
			G[k-1] = G[2+i];
			// Imaginary part:
			G[k] = - G[2+i+1]; 
		}
		return G; 
	}	 
	
	
	/*
	 *  for compatibility with vanilla ImageJ
	 *  twiddle factors are pre-computed 
	 */
	public static Pair<float[], float[]> fftR2Cp1d(float[] real, final int sign, final int nfft) {
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
		}
 		if (nfft<real.length) {
			real=FFTUtil.truncate(real,nfft);
		}  		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 		
 		float[] imag=new float[real.length]; 
		// Bit reversal:		
		FFTUtil.bit_reverse(real);
		final int len = real.length;
		
		//compute twiddle coefficients
 		Pair<double[], double[]> ptab=FFTUtil.expTable2 (len, sign);
 		final double[] cosTable =ptab.first;
 		final double[] sinTable =ptab.second;
			
		// Cooley-Tukey decimation-in-time radix-2 FFT
		int tablestep=len;
		
		for (int N = 2, hN=1; N <= len; N <<=1) {
			hN = N >>1;
			tablestep = tablestep >>1; //len >>c;
			//System.out.println("tablestep "+tablestep);
			for (int i = 0; i < len; i += N) {
				for (int j=i,  k = 0; j < i+hN; j++, k+= tablestep) {
					double tpre  =  real[j+hN] * cosTable[k] - imag[j+hN] * sinTable[k];
					double tpim  = real[j+hN] * sinTable[k] + imag[j+hN] * cosTable[k];
					real[j + hN] = (float)(real[j] - tpre);
					imag[j + hN] = (float)(imag[j] - tpim);
					real[j]+= tpre;
					imag[j]+= tpim;
				}
			}
		}
		
		if (sign==1) {
			for (int i=0; i<len; i++) {
				real[i] /=(float)len;
				imag[i] /=(float)len;
				if (abs(imag[i]) <=tol)
					imag[i]=0;
			}
		}
		return new Pair<float[], float[]>(real,imag);
	} // end
   
	
	/*
	 *  for compatibility with vanilla ImageJ
	 *  twiddle factors are precomputed 
	 */
	public static Pair<float[], float[]> fftR2Cq1d(float[] real, final int sign, final int nfft) {
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
		}
 		if (nfft<real.length) {
			real=FFTUtil.truncate(real,nfft);
		} 		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 		
 		float[] imag=new float[real.length];
 
		// Bit reversal:
		final int len = real.length;

		FFTUtil.bit_reverse(real);
		
		//compute twiddle coefficients
 		Pair<double[], double[]> ptab=FFTUtil.expTable (len, sign);
 		final double[]cosTable=ptab.first;
 		final double[]sinTable =ptab.second;
 		
		//compute twiddle coefficients
		// Cooley-Tukey decimation-in-time radix-2 FFT
		int c=1;
		for (int N = 2, hN=1; N <= len; N <<=1) {
			hN = N >>1;
			int tablestep = len >>c;
			for (int i = 0; i < len; i += N) {
				for (int j = i, k = 0; j < i + hN; j++, k += tablestep) {
					double tpre  =  real[j+hN] * cosTable[k] + imag[j+hN] * sinTable[k];
					double tpim  = -real[j+hN] * sinTable[k] + imag[j+hN] * cosTable[k];
					real[j + hN] = (float) (real[j] - tpre);
					imag[j + hN] = (float) (imag[j] - tpim);
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			c++;
		}

		if (sign==1) {
			for (int i=0; i<len; i++) {
				real[i] /=(float)len;
				imag[i] /=(float)len;
				if (abs(imag[i]) <=tol)
					imag[i]=0;
			}
		}
		
		return new Pair<float[], float[]>(real,imag);
	} // end
   

 
	/*
	 * computes forward FFT, real only input, size PoW2
	 */
	public static float[] rfft(float[] g){		
		int nfft=DSP.nfft(g.length); 
		// checks if the length of the array is power of two
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes interlaced FFT of the input
	    cfft(g, true);
	    //Util.printFloatArray(g);
	    // computes the final twiddle coefficients
	    double[] ct = FFTUtil.dPrecompute2(nfft);
	    // recalculates the values to get the rFFT (only the left part)
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]   
	    float[] GL = twiddleAndUnfold(ct,g);
		return GL;
	}
	
	public static double[] rfft(double[] g){		
		int nfft=nfft(g.length); 
		// checks if the length of the array is power of two
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes the final twiddle coefficients
	    double[] ct=dPrecompute2(nfft);
	    // computes interlaced FFT of the input
	    cfft(g, true);
	    //Util.printFloatArray(g);
	    // recalculates the values to get the rFFT (only the left part)
	    double[] GL = twiddleAndUnfold(ct,g);
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
		return GL;
	}
	
 
	public static float[] rfftp (float[] g,Pair<double[],double[]> ptab){
		
		// checks if the length of the array is power of two
		int nfft=DSP.nfft(g.length); 
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	 
	    // compute the FFT of the input (that's a black box for the method)
	    cfftp(g, true, ptab);
	   //System.out.println("cfftp step");
	   //Util.printFloatArray(g);
	   
	   	// computes the coefficients
	    double[] ct = FFTUtil.dPrecompute2(nfft);    
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
	    float[]  GL = twiddleAndUnfold(ct,g);
		return GL;
	}

	public static double[] rfftp (double[] g, Pair<double[],double[]> ptab){
		
		// checks if the length of the array is power of two
		int nfft=DSP.nfft(g.length); 
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	 
	    // compute the FFT of the input (that's a black box for the method)
	    cfftp(g, true, ptab);
	   //System.out.println("cfftp step");
	   //Util.printFloatArray(g);
	   
	   	// computes the coefficients
	    double[] ct = FFTUtil.dPrecompute2(nfft);    
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
	    double[]  GL = twiddleAndUnfold(ct,g);
		return GL;
	}
	
	
/*
 * Danielson and Lanczos FFT algorithm
 *  x must be interleaved real/complex array of power of 2 length
 * 
 */
 	public static void cfft( float[] x, boolean forward ) {
 		double theta, wpr, wpi, wr, wi, tmp ;	   
		int sgn=1, j;
		double rtemp, itemp ;
		
		int nfft=DSP.nfft(x.length);	
		if (x.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + x.length);

		int  len = x.length ; // ND
		if (forward) sgn=-1;
        bit_reverse2(x) ;
        
	    for(int N = 2, hN=1 ; N < len ;	N = hN ){
	         hN = N<<1 ;
	    	//System.out.println("N// "+N+" hn "+ hN);
	        theta =  sgn*TWOPI/N ;
	        tmp =  sin(0.5*theta );
	        wpr =  -2.0*tmp*tmp ;
	      //  wpr = cos(theta)-1.0;
	        wpi = sin(theta)  ;
	        wr = 1.0 ;
	        wi = 0.0 ;       	
	        for(int m = 0 ; m < N ; m += 2 ){
	        	//System.out.println("m "+ m);	        	
	            for(int i = m; i < len; i+= hN ){
	                j = i + N ;
	                //System.out.println("i "+ i+ " j "+j+" hn "+ hN+ " N "+ N);
	                /* complex multiplication */
 	                rtemp =  wr*x[j] - wi*x[j+1] ;
 	                itemp =  wr*x[j+1] + wi*x[j] ;	 
	                x[j] = (float) (x[i] - rtemp) ;
	                x[j+1] = (float) (x[i+1] - itemp) ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;            
	            }
	            /* complex multiplication */
	            tmp = wr;
				wr += tmp*wpr - wi*wpi ;
	            wi += wi*wpr + tmp*wpi ;
	        }
	    }
	 
	    // scale output
	    if (!forward) {
			for (int u=0; u<len; u++) {
				x[u] =(float) (2.0*x[u]/len);		
				if (Math.abs(x[u]) <=tol)
					x[u]=0;
			}
		}
	    
	}

public static void cfft( double[] x, boolean forward ) {
		
		int nfft=nfft(x.length);	
		if (x.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + x.length);

		int  len = x.length ; // ND
		double theta, wpr, wpi, wr, wi, tmp ;	   
		int sgn=1, j;
		if (forward) 
			sgn=-1;

        bit_reverse2(x) ;
        //System.out.println("len "+len);
 
        float rtemp, itemp ;
	    for(int N = 2, hN=1 ; N < len ;	N = hN ){
	    
	         hN = N<<1 ;
	    	//System.out.println("N// "+N+" hn "+ hN);
	    	
	        theta =  sgn*TWOPI/N ;
	        tmp= sin(0.5*theta );
	        wpr =  -2.*tmp*tmp ;
	        wpi =  sin(theta) ;
	        wr = 1.f ;
	        wi = 0.f ;       	
	        for(int m = 0 ; m < N ; m += 2 ){
	        	//System.out.println("m "+ m);	        	
	            for(int i = m; i < len; i+= hN ){
	                j = i + N ;
	                //System.out.println("i "+ i+ " j "+j+" hn "+ hN+ " N "+ N);
 	                rtemp = (float) (wr*x[j] - wi*x[j+1]) ;
 	                itemp = (float) (wr*x[j+1] + wi*x[j]) ;	 
	                x[j] = x[i] - rtemp ;
	                x[j+1] = x[i+1] - itemp ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;            
	            }
	            tmp = wr;
				wr += tmp*wpr - wi*wpi ;
	            wi += wi*wpr + tmp*wpi ;
	        }
	    }
	 
	    // scale output
	    if (!forward) {
			for (int u=0; u<len; u++) {
				x[u] /=(float)len;		
				if (Math.abs(x[u]) <=tol)
					x[u]=0;
			}
		}
	    
	}

	public static Pair<double[],double[]> cfftp(float[] x, boolean forward, Pair<double[],double[]> ptab){
	   
		int nfft=DSP.nfft(x.length);	
	    if (x.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + x.length);
	    
	    nfft=nfft>>1;
	    int len = x.length ;

	  //  System.out.println("len "+ len);
	    FFTUtil.bit_reverse2( x ) ;
	    int sign=forward? -1 : 1 ;
	    
		//compute twiddle coefficients
	    if (ptab==null)
	    	ptab=FFTUtil.expTable2 (len, sign);
 		 double[]cosTable=ptab.first;
 		 double[]sinTable=ptab.second;
 		
 		 if (cosTable.length!=nfft) {
 			 System.out.println("table length  mismatch" + cosTable.length );
 			ptab=FFTUtil.expTable (len, sign);
 	 		cosTable=ptab.first;
 	 		sinTable =ptab.second;
 		 }
 		
/* 		System.out.println("tab len " + cosTable.length);
 		ComplexArray icm=new ComplexArray(cosTable, sinTable, false);
        System.out.println(icm);*/
 		
 		int delta =0;  
 		int tablestep = x.length >>1;
	    for(int mmax = 2 ; mmax < len ; mmax = delta ){
	        delta = mmax<<1 ;
       // System.out.println ( " delta: " +delta +" tablestep "+ tablestep +" mmax "+ mmax);
	        int k=-tablestep;
	        for(int m = 0 ; m < mmax ; m += 2 ){
	        	k+=tablestep;
	        	double rtemp;
				//System.out.println (m+ " wr: " +wr +" wi "+ wi +" k "+k);
	            double itemp ;
	            for(int i = m; i < len ; i += delta ){
	                int j = i + mmax ;
	                // complex multiplication twiddle
	                rtemp =   (cosTable[k]*x[j] + sinTable[k]*x[j+1]) ;
	                itemp =   (cosTable[k]*x[j+1] - sinTable[k]*x[j]) ;
	                x[j] = (float) (x[i] - rtemp) ;
	                x[j+1] = (float) (x[i+1] - itemp) ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;            
	            }
	            
	        }
	        tablestep=tablestep>>1;
	    }
	 

	    if (!forward) {
			for (int u=0; u<len; u++) {
				x[u] =(float) (2.0*x[u]/len);			
				if (Math.abs(x[u]) <=tol)
					x[u]=0;
			}
		}
	    return ptab;
	}
	
	public static Pair<double[],double[]> cfftp(double[] x, boolean forward, Pair<double[],double[]> ptab){
		   
		int nfft=DSP.nfft(x.length);	
	    if (x.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + x.length);
	    
	    nfft=nfft>>1;
	    int len = x.length ;

	  //  System.out.println("len "+ len);
	    FFTUtil.bit_reverse2( x ) ;
	    int sign=forward? -1 : 1 ;
	    
		//compute twiddle coefficients
	    if (ptab==null)
	    	ptab=FFTUtil.expTable2 (len, sign);
 		 double[]cosTable=ptab.first;
 		 double[]sinTable=ptab.second;
 		
 		 if (cosTable.length!=nfft) {
 			 System.out.println("table length  mismatch" + cosTable.length );
 			ptab=FFTUtil.expTable (len, sign);
 	 		cosTable=ptab.first;
 	 		sinTable =ptab.second;
 		 }
 		
/* 		System.out.println("tab len " + cosTable.length);
 		ComplexArray icm=new ComplexArray(cosTable, sinTable, false);
        System.out.println(icm);*/
 		
 		int delta =0;  
 		int tablestep = x.length >>1;
	    for(int mmax = 2 ; mmax < len ; mmax = delta ){
	        delta = mmax<<1 ;
       // System.out.println ( " delta: " +delta +" tablestep "+ tablestep +" mmax "+ mmax);
	        int k=-tablestep;
	        for(int m = 0 ; m < mmax ; m += 2 ){
	        	k+=tablestep;
	        	double rtemp;
				//System.out.println (m+ " wr: " +wr +" wi "+ wi +" k "+k);
	            double itemp ;
	            for(int i = m; i < len ; i += delta ){
	                int j = i + mmax ;
	                // complex multiplication twiddle
	                rtemp =   (cosTable[k]*x[j] + sinTable[k]*x[j+1]) ;
	                itemp =   (cosTable[k]*x[j+1] - sinTable[k]*x[j]) ;
	                x[j] = (float) (x[i] - rtemp) ;
	                x[j+1] = (float) (x[i+1] - itemp) ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;            
	            }
	            
	        }
	        tablestep=tablestep>>1;
	    }
	 

	    if (!forward) {
			for (int u=0; u<len; u++) {
				x[u] =(float) (2.0*x[u]/len);			
				if (Math.abs(x[u]) <=tol)
					x[u]=0;
			}
		}
	    return ptab;
	}

	
	public static void fftC2C4(float[] X, float[] Y) {
		int N=X.length;
	    // N = 4 ^ M
	    int N1,N2;
	    int I1, I2, I3;
	    double CO1,CO2,CO3,SI1,SI2,SI3, A,B,C,theta;
	    double R1,R2,R3,R4, S1,S2,S3,S4;

	    N2 = N; I2 = 0; I3 = 0;
	    for (int k=0; k<N; k++) {
	        N1 = N2;
	        N2 = N2 / 4;
	        theta = TWOPI /N1;
	        A = 0.0;
	        for (int j=0; j < N2; j++) {
	            A = j*E;
	            
	            //Should be pre-calculated for optimization
	            CO1 = cos(A);
	            SI1 = sin(A);
	            
	            //B = A + A; twiddle a complex multiplication of a double argument
	            //CO2 = cos(B);  // CO2 = cos(A)*cos(A) -sin(A)*sin(A);
	            CO2= CO1*CO1 - SI1*SI1;
	            //SI2 = sin(B);  //SI2=2*cos(A)*sin(A);
	            SI2=2*CO1*SI1;
	            //  C = A + B; twiddle a complex multiplication
	            //CO3= cos(A)*cos(B)-sin(A)*sin(B);
	            CO3= CO1*CO2 - SI1*SI2;
	            //SI3 = cos(A)*sin(B)+sin(A)*cos(B)
	            SI3= CO1*SI2 + CO2*SI1;
	            
	            for (int i = j; i<N; i+=N1) {
	                I1 = i + N2;
	                I2 = I1 + N2;
	                I3 = I2 + N2;
	                R1 = X[i] + X[I2];
	                R3 = X[i] - X[I2];
	                S1 = Y[i] + Y[I2];
	                S3 = Y[i] - Y[I2];
	                R2 = X[I1] + X[I3];
	                R4 = X[I1] - X[I3];
	                S2 = Y[I1] + Y[I3];
	                S4 = Y[I1] - Y[I3];
	                X[i] = (float) (R1 + R2);
	                R2 = R1 - R2;
	                R1 = R3 - S4;
	                R3 = R3 + S4;
	                Y[i] = (float) (S1 + S2);
	                S2 = S1 - S2;
	                S1 = S3 + R4;
	                S3 = S3 - R4;
	                X[I1] = (float) (CO1*R3 + SI1*S3);
	                Y[I1] = (float) (CO1*S3 - SI1*R3);
	                
	                X[I2] = (float) (CO2*R2 + SI2*S2);
	                Y[I2] = (float) (CO2*S2 - SI2*R2);
	                
	                X[I3] = (float) (CO3*R1 + SI3*S1);
	                Y[I3] = (float) (CO3*S1 - SI3*R1);
	            } // end for i
	            //A+=theta;
	        } // end for j
	    } // end for k

	    // Radix-4 bit-reverse
	    float T;
	    int j = 0;
	    N2 = N>>2;
	    for (int i=0; i < N-1; i++) {
	        if (i < j) {
	            T = X[i];
	            X[i] = X[j];
	            X[j] = T;
	            T = Y[i];
	            Y[i] = Y[j];
	            Y[j] = T;
	        }
	        N1 = N2;
	        while ( j >= 3*N1 ) {
	            j -= 3*N1;
	            N1 >>= 2;
	        }
	        j += N1;
	    }
	} //////////
	
	 
	
	/* Maxima code
	 fft_rec(x, w):=block ( [n:length(x), odd, even, r, l, i, hlen, v, tmp ],
	local(v),
	if n=1 then 
		return (x),
	if oddp(n) then error ("length must be power of 2; found ", n),
	r:makelist(x[i], i, 1, length(x), 2),
	even: fft_rec (r, w^2),
	l:makelist(x[i],  i, 2, length(x), 2),
	odd: fft_rec (l, w^2),
	hlen: fix(n/2),
	for i:1 thru hlen do (
		tmp:demoivre(w^i)* odd[i],
		v[i]:expand(even[i] + tmp),
		v[i+hlen]:expand(even[i] - tmp)
	),
	listarray (v)	
	)$
	 */
	
	
	
	public static float[] zeroPaddEnd(float s[]) {
		final int n=s.length;
		final int k=DSP.nextpow2(n);
		final int n2=DSP.pow2(k);
		//System.out.println("n2 " +n2);
		float[] spad=new float[n2];
		System.arraycopy(s, 0, spad, 0, n);
		return spad;
	}
	
	public static double[] zeroPaddEnd(double s[]) {
		final int n=s.length;
		final int k=DSP.nextpow2(n);
		final int n2=DSP.pow2(k);
		//System.out.println("n2 " +n2);
		double[] spad=new double[n2];
		System.arraycopy(s, 0, spad, 0, n);
		return spad;
	}
	
	public static double[] zeroPaddMiddle(double s[]) {
		final int n=s.length;
		final int k=DSP.nextpow2(n);
		final int n2=DSP.pow2(k);
		//System.out.println("n2 " +n2);
		double[] spad=new double[n2];
		System.arraycopy(s, (n2-k)/2, spad, 0, n);
		return spad;
	}
	
	public static float[] zeroPaddMiddle(float s[]) {
		final int n=s.length;
		final int k=DSP.nextpow2(n);
		final int n2=DSP.pow2(k);
		//System.out.println("n2 " +n2);
		float[] spad=new float[n2];
		System.arraycopy(s, (n2-k)/2, spad, 0, n);
		return spad;
	}
	
	 /*
	  *  from PixLib 2014
	  */
	
	
	/**
	 * 
	 * @param arr
	 */
	public static void fftshift1d(double[] arr) {
		final int k=arr.length/2;		
		for (int i=0; i<k; i++) {
			final double tmp=arr[i];
			arr[i]=arr[k+i];
			arr[k+i]=tmp;			
		}
	}
	
	/**
	 * assumes IComplexArray structure
	 * @param arr
	 */
	public static void fftshift1c(double[] arr) {
		final int k=arr.length/2;		
		for (int i=0; i<k; i+=2) {
			final double tmpr=arr[i];
			final double tmpc=arr[i+1];
			arr[i]=arr[k+i];
			arr[i+1]=arr[k+i+1];
			arr[k+i]=tmpr;
			arr[k+i+1]=tmpc;
		}
	}
	
	public static void ifftshift1d(double[] arr) {
		final int k= arr.length/2;
		 
		for (int i= arr.length-1; i>=k; i--) {
			final double tmp=arr[i];
			arr[i]=arr[i-k];
			arr[i-k]=tmp;
		}
	}
	
	public static void ifftshift1c(double[] arr) {
		final int k= arr.length/2;
		 
		for (int i= arr.length-2; i>=k; i-=2) {
			final double tmpr=arr[i];
			final double tmpc=arr[i+1];
			arr[i]=arr[i-k];
			arr[i+1]=arr[i-k+1];
			arr[i-k]=tmpr;
			arr[i-k+1]=tmpc;
		}
	}	

	
	/**
	 * 
	 * @param arr
	 */
	public static void fftshift1d(float[] arr) {
		final int k=arr.length/2;		
		for (int i=0; i<k; i++) {
			final float tmp=arr[i];
			arr[i]=arr[k+i];
			arr[k+i]=tmp;			
		}
	}
	
	/**
	 * assumes IComplexArray structure
	 * @param arr
	 */
	public static void fftshift1c(float[] arr) {
		final int k=arr.length/2;		
		for (int i=0; i<k; i+=2) {
			final float tmpr=arr[i];
			final float tmpc=arr[i+1];
			arr[i]=arr[k+i];
			arr[i+1]=arr[k+i+1];
			arr[k+i]=tmpr;
			arr[k+i+1]=tmpc;
		}
	}
	
	public static void ifftshift1d(float[] arr) {
		final int k= arr.length/2;
		 
		for (int i= arr.length-1; i>=k; i--) {
			final float tmp=arr[i];
			arr[i]=arr[i-k];
			arr[i-k]=tmp;
		}
	}
	
	public static void ifftshift1c(float[] arr) {
		final int k= arr.length/2;
		 
		for (int i= arr.length-2; i>=k; i-=2) {
			final float tmpr=arr[i];
			final float tmpc=arr[i+1];
			arr[i]=arr[i-k];
			arr[i+1]=arr[i-k+1];
			arr[i-k]=tmpr;
			arr[i-k+1]=tmpc;
		}
	}	
	

	
} // end of class
///////////////////////////////
