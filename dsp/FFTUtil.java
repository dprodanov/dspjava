package dsp;


import static java.lang.Math.cos;
import static java.lang.Math.sin;
import ijaux.datatype.Pair;
import static dsp.FFTProc.TWOPI;

public class FFTUtil {

	public FFTUtil() {
		// TODO Auto-generated constructor stub
	}

	/*
	 * computes twiddle coefficients as a complex vector
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

	/*
	 * computes twiddle coefficients as a complex vector
	 *  quadrant symmetry
	 */
	public static Pair<double[], double[]> expTable2(int len, int sign) {
		int hlen = len>>1;
		int qlen = len>>2;
		final double[]cosTable = new double[hlen];
		final double[]sinTable = new double[hlen];
		
		final double theta =  TWOPI/(double)len;
		
		cosTable[0]=1;
		// we exploit the symmetries of the 
		// complex roots of unity
		if (sign <0) {
			for (int i = 1; i<qlen; i++) {
				double c=cos(i*theta);
				cosTable[i     ] = c;
				sinTable[qlen-i] = c;
				
				cosTable[hlen-i] = -c;
				sinTable[qlen+i] = c;
			}
			sinTable[qlen]=cosTable[0];
		} else {
			for (int i = 1; i<qlen; i++) {
				double c=cos(i*theta);
				cosTable[i     ] = c;
				sinTable[qlen-i] = -c;
				
				cosTable[hlen-i] = -c;
				sinTable[qlen+i] = -c;
			}
			sinTable[qlen]=-cosTable[0];
		}
		
		
		return Pair.of(cosTable,sinTable);
	}

	/*
	 * computes twiddle coefficients interlaced
	 */
	public static double[] iexpTable (int len, int sign) {
		final double[] wtable = new double[len];	
		final double theta = -sign*TWOPI/(double)len;
		wtable[0]=1;
		double arg=theta;		
		for (int i = 2; i < len; i+=2) {
			wtable[i] = cos(arg);
			wtable[i+1] = sin(arg );
			arg += theta ;
		}
		return wtable;
	}

	/*
	 *  computes twiddle coefficients as a complex vector
	 *  abscissa symmetry; coefficients interlaced
	 */
	public static double[] iexpTable2 (int len, int sign) {
		
		final double[] wtable = new double[len];	
		final double theta = TWOPI/(double)len;
		
		int qlen = len>>1;
		wtable[0]=1;
		if (sign <0 ) {
			for (int i = 2, j=1; i<qlen; i+=2, j++) {
				double c=cos(j*theta);
					   wtable[i] = c;
				wtable[qlen-i+1] = c;
				wtable[len-i]   = -c;
				wtable[qlen+i+1] = c;
			}
			wtable[qlen+1]=wtable[0];
		} else {
			for (int i = 2, j=1; i<qlen; i+=2, j++) {
				double c=cos(j*theta);
					   wtable[i] = c;
				wtable[qlen-i+1] = -c;
				wtable[len-i]   = -c;
				wtable[qlen+i+1] = -c;
			}
			wtable[qlen+1]=-wtable[0];
		}
		return wtable;
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

	/*
	 *   bitreverse places float array x containing N/2 complex values
	 *   into bit reversed indexing
	 */
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
		 * @param data
		 */
		public static void printvector(float[] data) {
			for (int i=0; i<data.length; i++) {		
					System.out.print(data[i]+",");		
			}
		}
		
		public static void printvector(double[] data) {
			for (int i=0; i<data.length; i++) {		
					System.out.print(data[i]+",");		
			}
		}

		public static float[] truncate(float s[], int k) {
			if (k>s.length)
				throw new IllegalArgumentException(k+">"+s.length);
			float[] spad=new float[k];
			System.arraycopy(s, 0, spad, 0, k);
			return spad;
		}

		/* 
		 * Method for computing the sine and cosine coefficients,
		 * Even indicies contain the real part and odd indicies contain the imaginary part
		 * input:
		 *  N - Length of the array which contains real numbers, must be power of two.
		 * output:
		 *  A - matrix of the cos coefficients
		 */
		public static  double[]  dPrecompute2 (int N){
			double[] A = new double[N]; 
			//System.out.println(N);
			for (int i=0; i<N; i+=2){
				double theta=i*Math.PI/N;
				// Real part:
				A[i  ]=  0.5*(1.0 - sin(theta));
				// Imag part:
				A[i+1]= -0.5*cos(theta);
			}
			return A;
		}

}
