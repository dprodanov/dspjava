package dsp;


import ij.process.*;
import ijaux.Util;
import ijaux.datatype.*;
import ijaux.dsp.DSP;
import ijaux.dsp.FFT;
import ijaux.scale.IJLineIteratorIP;
import static java.lang.Math.*;
import static ijaux.dsp.FFT.TWOPI;

/*
 *  FFT radix 2 implementations
 *  (C) Dimiter Prodanov, 
 *  (C) Tomasz Tomasz Konopczy\'nski
 */
public class FFTProc3 {
	
	private static final float tol=2E-15f;
	
	double[] cosTable=null, sinTable=null, ctable=null;
	
	private int len=-1;
	
	public FFTProc3(int n) {
		n=FFT.nfft(n);
		init(n, -1);
		len=n;		
	}

	/**
	 * @param n
	 */
	private void init(int n,  int sign) {
		Pair<double[], double[]> ptable=expTable(n, sign);
		cosTable=ptable.first;
		sinTable=ptable.second;
		ctable=dPrecompute2(n);
	}
 
	public Pair<FloatProcessor,FloatProcessor> fftC2C1D(Pair<FloatProcessor,FloatProcessor> pr, int sign, int xdir) {
			
		FloatProcessor rp=pr.first;
		FloatProcessor cp=pr.second;
		IJLineIteratorIP<float[]> iter_r= new IJLineIteratorIP<float[]>(rp, xdir);
		IJLineIteratorIP<float[]> iter_c= new IJLineIteratorIP<float[]>(cp, xdir);
		
		int width=rp.getWidth();
		int height=rp.getHeight();
		
		int[] dim={width, height};
		
		int bs=dim[xdir];
		
		//Util.printIntArray(dim);
		int nfft=FFT.nfft(bs);
		
		//System.out.println ("\nbs " +bs +" "+DSP.nextpow2(bs)+" nfft "+nfft);
		if (xdir==0)
			width=nfft;
		if (xdir==1)
			height=nfft;
	 	
		//System.out.println ("width " + width +" height "+ height);
		
		FloatProcessor re=new FloatProcessor(width, height);
		FloatProcessor im=new FloatProcessor(width, height);
		
		IJLineIteratorIP<float[]> iter2= new IJLineIteratorIP<float[]>(new int[]{width, height}, 16, xdir);
		
		int cnt=0;
 
		while (iter_r.hasNext() && iter_c.hasNext() ) {
			//System.out.println(" c: "+cnt);
			final float[] real=iter_r.next();	
			final float[] imag=iter_c.next();	
 
		 	final Pair<float[], float[]> carr= fftC2C1d(real, imag, sign, nfft);
			
		/*	Util.printFloatArray(carr.first);
			Util.printFloatArray(carr.second);			
			System.out.println(); */
			
			iter2.putLineFloat(re,carr.first, cnt, xdir);
			iter2.putLineFloat(im,carr.second, cnt, xdir);
			cnt++;		
		}
	
		return Pair.of(re,im);
	}
	
	public Pair<FloatProcessor,FloatProcessor> fftR2C1D(FloatProcessor fp, int sign, int xdir) {
		IJLineIteratorIP<float[]> iter= new IJLineIteratorIP<float[]>(fp, xdir);
		
		int width=fp.getWidth();
		int height=fp.getHeight();		
		int[] dim={width, height};		
		int bs=dim[xdir];
		
		//Util.printIntArray(dim);
		int nfft=FFT.nfft(bs);
		
		//System.out.println ("\nbs " +bs +" "+DSP.nextpow2(bs)+" nfft "+nfft);
		if (xdir==0)
			width=nfft;
		if (xdir==1)
			height=nfft;
	 	
		//System.out.println ("width " + width +" height "+ height);
		
		FloatProcessor re=new FloatProcessor(width, height);
		FloatProcessor im=new FloatProcessor(width, height);
		
		IJLineIteratorIP<float[]> iter2= new IJLineIteratorIP<float[]>(new int[]{width, height}, 16, xdir);
		int cnt=0;
		while (iter.hasNext()) {
			//System.out.println(" c: "+cnt);
			final float[] real=iter.next();	
			 
		 	final Pair<float[], float[]> carr= fftR2C1d(real, sign, nfft);
/*			
		 	Util.printFloatArray(carr.first);
			Util.printFloatArray(carr.second);			
			System.out.println(); */
			
			iter2.putLineFloat(re,carr.first, cnt, xdir);
			iter2.putLineFloat(im,carr.second, cnt, xdir);
			cnt++;		
		}
	
		return Pair.of(re,im);
	}
	
	private final int Ox=0, Oy=1, Oz=2, Oc=3, Ot=4;
	
	public static boolean debug=false;
	
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
	
	public Pair<FloatProcessor,FloatProcessor> fftR2C2D(FloatProcessor fp, int sign) {
 
		int[] dim=nfftp(fp);
		int nfftx=dim[0];
		int nffty=dim[1];
 		
		//System.out.println ("width " + nfftx +" height "+ nffty);
		
		FloatProcessor re=new FloatProcessor(nfftx, nffty);
		FloatProcessor im=new FloatProcessor(nfftx, nffty);
		
		IJLineIteratorIP<float[]> iterin_x= new IJLineIteratorIP<float[]>(fp, Ox);
		
		//IJLineIteratorIP<float[]> iterout_x= new IJLineIteratorIP<float[]>(new int[]{width, height},32 ,Ox);
		IJLineIteratorIP<float[]> iterout_x= new IJLineIteratorIP<float[]>(re, Ox);
		
		
		int cnt=0;
 
		while (iterin_x.hasNext()) {
			//System.out.println(" c: "+cnt);
			float[] real=iterin_x.next();	 
		 	final Pair<float[], float[]> carr= fftR2C1d(real, sign, nfftx);
			
/*		 	Util.printFloatArray(carr.first);
			Util.printFloatArray(carr.second);			
			System.out.println(); 
			*/
			iterout_x.putLineFloat(re,carr.first, cnt, Ox);
			iterout_x.putLineFloat(im,carr.second, cnt, Ox);
			cnt++;		
		}
		
		
		//IJLineIteratorIP<float[]> iterin_xre= new IJLineIteratorIP<float[]>(new int[]{width, height}, 32, Oy);
		IJLineIteratorIP<float[]> iterin_xre= new IJLineIteratorIP<float[]>(re, Oy);
		IJLineIteratorIP<float[]> iterin_xim= new IJLineIteratorIP<float[]>(im, Oy);
		
/*		FloatProcessor re2=new FloatProcessor(width, height);
		FloatProcessor im2=new FloatProcessor(width, height)*/;
		//IJLineIteratorIP<float[]> iterout_y= new IJLineIteratorIP<float[]>(re, Oy);
		//IJLineIteratorIP<float[]> iterout_y2= new IJLineIteratorIP<float[]>(re, Oy);
		IJLineIteratorIP<float[]> iterout_y= new IJLineIteratorIP<float[]>(new int[]{nfftx, nffty}, 32, Oy);
		//iterout_x=new IJLineIteratorIP<float[]>(new int[]{width, height},16 ,Oy);
		cnt=0;
		while (iterin_xre.hasNext()) {
			//System.out.println(" c: "+cnt);
			//float[] real=iterin_xre.getLineFloat(re, cnt, Oy);
			float[] real=iterin_xre.next();
			float[] imag=iterin_xim.next();
			//float[] imag=iterin_xre.getLineFloat(im, cnt, Oy);
		 	final Pair<float[], float[]> carr= fftC2C1d(real, imag, sign, nffty);
			
		 	//Util.printFloatArray(carr.first);
			//Util.printFloatArray(carr.second);			
			//System.out.println(); 
			
			iterout_y.putLineFloat(re,carr.first, cnt, Oy);
			iterout_y.putLineFloat(im,carr.second, cnt, Oy);
			
			cnt++;
			
			//iterin_xre.fwd();
		}
		
		return Pair.of(re,im);
	}
	
	/*
	 *  for compatibility with vanilla ImageJ
	 */
	public static Pair<float[], float[]> fftC2C1d(float[] real, float[] imag, final int sign, final int nfft) {	
		
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
			imag=zeroPaddEnd(imag);
		}
 		if (nfft<real.length) {
			real=truncate(real,nfft);
			imag=truncate(imag,nfft);
		} 
		if (real.length!= imag.length)
			throw new IllegalArgumentException("Array size mismatch " + real.length +" " +imag.length );

		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 
		
		// Bit reversal:
		final int len = real.length;
		//final int hlen = len/2;
		bit_reverse(real);
		bit_reverse(imag);
		
		double theta, tmp, alpha, beta, wr, wi;
		float tmpr, tmpi;
		//System.out.println("len "+ real.length);
		// Danielson-Lanczos algorithm:
		for (int N = 2, hN = 1; N <= len; hN = N, N <<= 1) {
			//System.out.println("N // "+ N+ " hN "+hN);
			theta = sign*TWOPI/N;
			tmp = Math.sin(0.5*theta);
			alpha = -2.0*tmp*tmp;
			beta = Math.sin(theta);
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
			real=truncate(real,nfft);
		} 
 		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 
 		float[] imag=new float[real.length];
 
		final int len = real.length;
		bit_reverse(real);
		
		// Danielson-Lanczos algorithm:
		for (int N = 2, hN = 1; N <= len; hN = N, N <<= 1) {
			final double theta = sign*TWOPI/N;
			double tmp = sin(0.5*theta);
			final double alpha = -2.0*tmp*tmp;
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
	
	/*
	 * computes twiddle coefficients as a complex vector
	 *  abscissa symmetry
	 */
	public static Pair<double[], double[]> expTable(int len, int sign) {
		int hlen = len>>1;
		final double[]cosTable = new double[hlen];
		final double[]sinTable = new double[hlen];
		
		final double theta = - sign*TWOPI/(double)len;
		cosTable[0]=1;
		double arg=theta;		
		for (int i = 1; i<hlen; i++) {
			cosTable[i] = cos(arg);
			sinTable[i] = sin(arg );
			arg += theta ;
		}
		return Pair.of(cosTable,sinTable);
	}
	
	public static Pair<double[], double[]> expTable2(int len, int sign) {
		int hlen = len>>1;
		int qlen = len>>2;
		final double[]cosTable = new double[hlen];
		final double[]sinTable = new double[hlen];
		
		final double theta =  sign*TWOPI/(double)len;
		
		cosTable[0]=1;
		// we exploit the symmetries of the 
		// complex roots of unity
		for (int i = 1; i<qlen; i++) {
			double c=cos(i*theta);
			cosTable[i] = c;
			sinTable[qlen-i] = c;
			cosTable[hlen-i] = -c;
			sinTable[qlen+i] = c;
		}
		sinTable[hlen/2]=-cosTable[0];
		
		return Pair.of(cosTable,sinTable);
	}
	
	/*
	 * computes twiddle coefficients interlaced
	 */
	public static double[] iexpTable (int len, int sign) {

		final double[] wtable = new double[len];	
		final double theta = sign*TWOPI/(double)len;
		wtable[0]=1;
		double arg=theta;		
		for (int i = 2; i < len>>1; i+=2) {
			wtable[i] = cos(arg);
			wtable[i+1] = -sin(arg );
			arg += theta ;
		}
		return wtable;
	}
	
	/* 
	 * Method for computing the sine and cosine coefficients,
	 * Even indicies contain the real part and odd indicies contain the imaginary part
	 * input:
	 *  N - Length of the array which contains real numbers, must be power of two.
	 * output:
	 *  AB - matrix of the coefficients
	 */
	public static Pair<double[], double[]> dPrecompute (int N){
	
		double[] A = new double[N];
		double[] B = new double[N];
		//System.out.println(N);
		for (int i=0; i<N; i+=2){
			double theta=i*Math.PI/N;
			// Real part:
			A[i]=  0.5*(1.0 - sin(theta));
			B[i]=  1.0 - A[i];
			// Imag part:
			A[i+1]= - 0.5*cos(theta);
			B[i+1]= - A[i+1];
			//System.out.println(i+ " A " +A[i]+" "+A[i+1] );
			//System.out.println(theta+ " B " +B[i]+" "+B[i+1] );
		}
		return new Pair<double[], double[]>(A,B);
	}
	
	/* 
	 * Method for computing the sine and cosine coefficients,
	 * Even indicies contain the real part and odd indicies contain the imaginary part
	 * input:
	 *  N - Length of the array which contains real numbers, must be power of two.
	 * output:
	 *  AB - matrix of the coefficients
	 */
	public static  double[]  dPrecompute2 (int N){
		double[] A = new double[N]; 
		//System.out.println(N);
		for (int i=0; i<N; i+=2){
			double theta=i*Math.PI/N;
			// Real part:
			A[i]=  0.5*(1.0 - sin(theta));
			// Imag part:
			A[i+1]= - 0.5*cos(theta);
		}
		return A;
	}
	/* 
	 * Method for computing the sine and cosine coefficients,
	 * Even indices contain the real part and odd indices contain the imaginary part
	 * input:
	 *  N - Length of the array which contains real numbers, must be power of two.
	 * output:
	 *  AB - matrix of the coefficients
	 */
	public static Pair<float[], float[]> fPrecompute (int N){
		float[] A = new float[N];
		float[] B = new float[N];
		for (int i=0; i<(N/2); i++){
			double df=i*TWOPI/N;
			int i2=2*i;
			// Real part:
			A[i2]=  (float) (0.5*(1.0 - sin(df)));
			B[i2]=  (float) (0.5*(1.0 + sin(df)));
			// Imag part:
			A[i2+1]= - (float)(0.5*cos(df));
			B[i2+1]=   (float) (0.5*cos(df));			
		}
		return new Pair<float[], float[]>(A,B);
	}
	 
	/*
	 * Computes the true FFT values (G array) from the array X.
	 * input:
	 *  N - Length of the array which contains real numbers, must be power of two.
	 *  AB - matrix of the coeficients
	 *  X1 - output from the FFT. 
	 * output:
	 *  G - true FFT values.
	 */
	public static float[] ComputeGfromX (int N, Pair<float[], float[]>ABtable, float[] X){
		float[] G = new float[N+2];
		int hN = N/2;
		for (int k = 0; k < hN; k++){
			int k2= k*2;
			int u = N-k2;
			if (u+1>N) u=u-N;			
			//System.out.println((k2+1)+" "+ (u+1)+ " "+ (u));
			// Real part:
			G[k2  ]=(X[k2]*ABtable.first[k2] 
					- X[k2+1]*ABtable.first[k2+1]
					+ X[u]*ABtable.second[k2] + X[u+1]*ABtable.second[k2+1]);
			// Imag part:
			G[k2+1]=(X[k2+1]*ABtable.first[k2] 
					+ X[k2]*ABtable.first[k2+1]
					+ X[u]*ABtable.second[k2+1] - X[u+1]*ABtable.second[k2]);
		}		
		//last real
		G[N]=X[0]-X[1];
		//last imag
		G[N+1]=0;
		return G; 
	}
	
	/*
	 * Computes the true FFT values (G array) from the array X.
	 * input:
	 *  N - Length of the array which contains real numbers, must be power of two.
	 *  AB - matrix of the coefficients
	 *  X1 - output from the FFT. 
	 * output:
	 *  G - true FFT values.
	 */
	/*public static float[] ComputeGfromX2 (Pair<double[], double[]>ABtable, float[] X){
		int N=X.length;
		//System.out.println(N);
		float[] G = new float[N+2];
		for (int k = 0; k < N; k+=2){
			int u = N-k;
			//if (u+1>N) u=u-N;	
			if (1>k) u=u-N;	
			//System.out.println(k  +" "+u);
			//System.out.println(k+"  "+(k+1)+" "+ (u+1)+ " "+ (u)+" "+ (N-k-1));
			// Real part:
			G[k  ]=(float)(X[k]*ABtable.first[k] 
					- X[k+1]*ABtable.first[k+1]
					+ X[u]*ABtable.second[k] 
					+ X[u+1]*ABtable.second[k+1]);
			// Imaginary part:
			G[k+1]=(float)(X[k+1]*ABtable.first[k] 
					+ X[k]*ABtable.first[k+1]
					+ X[u]*ABtable.second[k+1] 
					- X[u+1]*ABtable.second[k]);
			//System.out.println("// "+G[k  ]+"  "+X[k] +" "+ X[k+1]+" "+X[u]+" "+X[u+1]);
			//System.out.println("// " +"  "+ABtable.first[k] +" "+ ABtable.first[k+1]+" "+ABtable.second[k]+" "+ABtable.second[k+1]);
			
			//System.out.println("// "+G[k  ]+"  "+G[k+1]);
		}		
		//last real
		G[N]=X[0]-X[1];
		//last imag
		G[N+1]=0;
		return G; 
	}*/
	
	/*public static float[] ComputeGfromX3 (Pair<double[], double[]>ABtable, float[] X){
		int N=X.length;
		//System.out.println(N);
		float[] G = new float[N<<1];
		for (int k = 0; k < N; k+=2){
			int u = N-k;
			//if (u+1>N) u=u-N;	
			if (1>k) u=u-N;	
			//System.out.println(k  +" "+u);
			//System.out.println(k+"  "+(k+1)+" "+ (u+1)+ " "+ (u)+" "+ (N-k-1));
			// Real part:
			G[k  ]=(float)(X[k]*ABtable.first[k] 
					- X[k+1]*ABtable.first[k+1]
					+ X[u]*ABtable.second[k] 
					+ X[u+1]*ABtable.second[k+1]);
			// Imaginary part:
			G[k+1]=(float)(X[k+1]*ABtable.first[k] 
					+ X[k]*ABtable.first[k+1]
					+ X[u]*ABtable.second[k+1] 
					- X[u+1]*ABtable.second[k]);
			//System.out.println("// "+G[k  ]+"  "+X[k] +" "+ X[k+1]+" "+X[u]+" "+X[u+1]);
			//System.out.println("// " +"  "+ABtable.first[k] +" "+ ABtable.first[k+1]+" "+ABtable.second[k]+" "+ABtable.second[k+1]);
			
			//System.out.println("// "+G[k  ]+"  "+G[k+1]);
		}		
		//last real
		G[N]=X[0]-X[1];
		//last imag
		G[N+1]=0;
		
		for (int i = 0; i<N; i+=2){
			int  k = 2*N -1-i;// full measure
			//System.out.println(i+" "+(k));
			// Real part:
			G[k-1] = G[2+i];
			// Imag part:
			G[k] = - G[2+i+1]; 
		}
		return G; 
	}*/
	
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
	/*
	 *  for compatibility with vanilla ImageJ
	 *  twiddle factors are precomputed 
	 */
	public static Pair<float[], float[]> fftR2Cp1d(float[] real, final int sign, final int nfft) {
		if (nfft>real.length) {
			real=zeroPaddEnd(real);
		}
 		if (nfft<real.length) {
			real=truncate(real,nfft);
		}  		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 		
 		float[] imag=new float[real.length]; 
		// Bit reversal:		
		bit_reverse(real);
		final int len = real.length;
		
		//compute twiddle coefficients
 		Pair<double[], double[]> ptab=expTable (len, sign);
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
			real=truncate(real,nfft);
		} 		
 		if (Math.abs(sign)!=1)
 			throw new IllegalArgumentException("illegal value "+sign);
 		
 		float[] imag=new float[real.length];
 
		// Bit reversal:
		final int len = real.length;
		//final int hlen = real.length/2;
		//System.out.println("hlen " +hlen);
		bit_reverse(real);
		
		//compute twiddle coefficients
 		Pair<double[], double[]> ptab=expTable (len, sign);
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
	 * Computes the rFFT from the array g.
	 * input:
	 *  g - array which contains real numbers only,
	 *  length of the array must be power of two.
	 * output:
	 *  G - rFFT of the array g.
	 */		
	public static float[] rfft(float[] g, boolean forward){
		int nfft=FFT.nfft(g.length); 
		// checks if the length of the array is power of two
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes the coefficients
	    Pair<float[], float[]> AB=fPrecompute(nfft); 
	    // computes the FFT of the input (that's a black box for the method)
	    cfft(g, forward);
	    // recalculates the values to get the rFFT (only the left part)
	    float[] GL = ComputeGfromX(nfft,AB,g);
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
	    float[] G  = unfoldFT(GL);
		return G;
	}
	
/*
 * computes forward FFT, real only input, size PoW2
 */
public static float[] rfft2(float[] g){		
		int nfft=FFT.nfft(g.length); 
		// checks if the length of the array is power of two
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes the final twiddle coefficients
	    double[] ct=dPrecompute2(nfft);
	    // computes interlaced FFT of the input
	    cfft(g, true);
	    //Util.printFloatArray(g);
	    // recalculates the values to get the rFFT (only the left part)
	    float[] GL = twiddleAndUnfold(ct,g);
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
		return GL;
	}
	
	/*
	 * Computes the rFFT from the array g.
	 * input:
	 *  g - array which contains real numbers only,
	 *  length of the array must be power of two.
	 * output:
	 *  G - rFFT of the array g.
	 */
	public static float[] rfftp (float[] g, boolean forward ,Pair<double[],double[]> ptab){
		
		// checks if the length of the array is power of two
		int nfft=FFT.nfft(g.length); 
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes the coefficients
	    Pair<float[], float[]> AB=fPrecompute(nfft); 
	    // compute the FFT of the input (that's a black box for the method)
	    cfftp(g, forward, ptab);
	    // recalculates the values to get the rFFT (only the left part)
	    float[] GL = ComputeGfromX(nfft,AB,g);
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
	    float[] G  = unfoldFT(GL);
		return G;
	}
	
public static float[] rfftp2 (float[] g,Pair<double[],double[]> ptab){
		
		// checks if the length of the array is power of two
		int nfft=FFT.nfft(g.length); 
	    if (g.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + g.length);
	    // computes the coefficients
	    double[] ct=dPrecompute2(nfft); 
	    // compute the FFT of the input (that's a black box for the method)
	   cfftp(g, true, ptab);
	   //System.out.println("cfftp step");
	   //Util.printFloatArray(g);

	    float[] GL = twiddleAndUnfold(ct,g);
	    // calculates the right part which is a mirror of the left [z(i)==z*(N-i)]
		return GL;
	}
	
/*
 * Danielson and Lanczos FFT algorithm
 *  x must be interleaved real/complex array of power of 2 length
 * 
 */
 	public static void cfft( float[] x, boolean forward ) {
		
		int nfft=FFT.nfft(x.length);	
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

	static float[] unfoldFT (float[] G){
		int N = G.length;
		float[] A = new float[2*N-4];
		System.arraycopy(G,0,A,0,N);
	
		int hm = N/2; // half measure
		
		// System.out.println(hm+" "+k);
		for (int i = 0; i<hm; i+=2){
			int  k = (2*N -5)-i;// full measure
			System.out.println(i+" "+(k));
			// Real part:
			A[k-1] = G[2+i];
			// Imag part:
			A[k] = - G[2+i+1]; 
		}
		return A;
	}
	

	public static Pair<double[],double[]> cfftp(float[] x, boolean forward, Pair<double[],double[]> ptab){
	   
		int nfft=FFT.nfft(x.length);	
	    if (x.length!= nfft)
			throw new IllegalArgumentException("Array size mismatch " + x.length);
	    
	    nfft=nfft>>1;
	    int len = x.length ;
	    //bit_reverse2( x, len ) ;
	  //  System.out.println("len "+ len);
	    bit_reverse2( x ) ;
	    int sign=forward? -1 : 1 ;
		//compute twiddle coefficients
	    if (ptab==null)
	    	ptab=expTable (len, sign);
 		 double[]cosTable=ptab.first;
 		 double[]sinTable=ptab.second;
 		
 		 if (cosTable.length!=nfft) {
 			ptab=expTable (len, sign);
 	 		cosTable=ptab.first;
 	 		sinTable =ptab.second;
 		 }
 		
/* 		System.out.println("tab len " + cosTable.length);
 		ComplexArray icm=new ComplexArray(cosTable, sinTable, false);
        System.out.println(icm);*/
 		
 		int  delta =0;  
 		int tablestep = x.length >>1;
	    for(int mmax = 2 ; mmax < len ; mmax = delta ){
	        delta = mmax<<1 ;
       // System.out.println ( " delta: " +delta +" tablestep "+ tablestep +" mmax "+ mmax);
	        int k=-tablestep;
	        for(int m = 0 ; m < mmax ; m += 2 ){
	        	k+=tablestep;
	        	//System.out.println (m+ " wr: " +wr +" wi "+ wi +" k "+k);
	            float rtemp, itemp ;
	            for(int i = m; i < len ; i += delta ){
	                int j = i + mmax ;
	                rtemp = (float) (cosTable[k]*x[j] + sinTable[k]*x[j+1]) ;
	                itemp = (float) (cosTable[k]*x[j+1] - sinTable[k]*x[j]) ;
	                x[j] = x[i] - rtemp ;
	                x[j+1] = x[i+1] - itemp ;
	                x[i] += rtemp ;
	                x[i+1] += itemp ;            
	            }
	            
	        }
	        tablestep=tablestep>>1;
	    }
	 

	    if (!forward) {
			for (int u=0; u<len; u++) {
				x[u] /=(float)len;		
				if (Math.abs(x[u]) <=tol)
					x[u]=0;
			}
		}
	    return ptab;
	}

	//-----------------------------------------------------------------------------
	// name: bit_reverse()
	// desc: bitreverse places float array x containing N/2 complex values
//	       into bit-reversed order
	//-----------------------------------------------------------------------------
	public static void bit_reverse2( float[] x) {	
		int N=x.length;
	    float rtemp, itemp ;
	    int i, j, m ;
	    for( i = j = 0 ; i < N ; i += 2, j += m )
	    {
	        if( j > i ){
	            rtemp = x[j] ; 
	            itemp = x[j+1] ; /* complex exchange */
	            x[j] = x[i] ; 
	            x[j+1] = x[i+1] ;
	            x[i] = rtemp ; 
	            x[i+1] = itemp ;
	        }

	        for( m = N>>1 ; m >= 2 && j >= m ; m >>= 1 )
	            j -= m ;
	    }
	}
	/**
	 * @param array
	 * @return
	 */
	public static void bit_reverse(float[] array) {
		final int len = array.length;
		final int hlen = len/2;
		//for (int i=0, j=0; i<len; ++i) {
		for (int i=0, j=0; i<len; i++) {
			if (j > i) {
				final float rt = array[j]; 
				array[j] = array[i]; 
				array[i] = rt;
			}
			int m = hlen;
			while (m >= 2 && j >= m) { 
				j -= m; 
				m >>= 1; 
	        }
			j += m;
		}
		//return len;
	}
	
	public static float[] arrD2F (double[] arr) {
		float[] ret=new float[arr.length];		
		for (int i=0; i<ret.length; i++)
			ret[i]=(float) arr[i];		
		return ret;
	}

	public static double[] arrF2D (float[] arr) {
		double[] ret=new double[arr.length];	
		for (int i=0; i<ret.length; i++)
			ret[i]=arr[i];		
		return ret;
	}
	
	public static float[] truncate(float s[], int k) {
		if (k>s.length)
			throw new IllegalArgumentException(k+">"+s.length);
		float[] spad=new float[k];
		System.arraycopy(s, 0, spad, 0, k);
		return spad;
	}
	
	public static float[] zeroPaddEnd(float s[]) {
		final int n=s.length;
		final int k=DSP.nextpow2(n);
		final int n2=DSP.pow2(k);
		//System.out.println("n2 " +n2);
		float[] spad=new float[n2];
		System.arraycopy(s, 0, spad, 0, n);
		return spad;
	}
	
	/**
	 * @param data
	 */
	private static void printvector(float[] data) {
		for (int i=0; i<data.length; i++) {		
				System.out.print(data[i]+",");		
		}
	}
	
	/*	
	 *  for compatibility with vanilla ImageJ
	 
	public static ComplexArray fft1d(float[] r, float[] im, final int sign, final int nfft) {
		//double[] imag=new double[real.length];
		//System.out.println("sz real: " + real.length);
		//System.out.println("sz imag: " + imag.length);
		double[] real=arrF2D(r);
		double[] imag=arrF2D(im);
		return FFT.fft1d(real, imag, sign, nfft);
	} // end
*/	
}
