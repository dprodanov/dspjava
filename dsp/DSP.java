package dsp;
import Jama.Matrix;


public class DSP {
	
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
	public static float[] arrF2Ci (float[] arr) {
		float[] ret=new float[2*arr.length];
		for (int i=0, k=0; i<arr.length; i++, k+=2) 
			ret[k]= arr[i];
		return ret;
	}
	
	public static double[] arrF2Cid (float[] arr) {
		double[] ret=new double[2*arr.length];
		for (int i=0, k=0; i<arr.length; i++, k+=2) 
			ret[k]= arr[i];
		return ret;
	}

	public static double[] arrD2Ci (double[] arr) {
		double[] ret=new double[2*arr.length];	
		for (int i=0, k=0; i<arr.length; i++, k+=2) 
			ret[k]= arr[i];
		return ret;
	}
	
	
	// implements Matlab's function flipr
		public static void flipr (double[][] kernel) {		
			for (int i=0; i< kernel.length; i++) {
				flip( kernel[i]);
			}
		}
		
		public static void flipud (Matrix m) {
			flipud (m.getArray());
		}
		
		// implements Matlab's function flipud
		public static void flipud (double[][] kernel) {		
			final int s=kernel.length-1;
			for (int i=0; i< kernel.length/2; i++) {
				final double[] c=kernel[i];
				kernel[i]=kernel[s-i];
				kernel[s-i]=c;
			}
			
		}
		
		// implements Matlab's function flipud
		public static void flipud (float[][] kernel) {		
			final int s=kernel.length-1;
			for (int i=0; i< kernel.length/2; i++) {
				final float[] c=kernel[i];
				kernel[i]=kernel[s-i];
				kernel[s-i]=c;
			}
			
		}
		/**
		 * @param kernel
		 */
		public static void flip(double[] kernel) {
			final int s=kernel.length-1;
			for (int i=0; i< kernel.length/2; i++) {
				final double c=kernel[i];
				kernel[i]=kernel[s-i];
				kernel[s-i]=c;
			}
		}
		
		public static void flip(float[] kernel) {
			final int s=kernel.length-1;
			for (int i=0; i< kernel.length/2; i++) {
				final float c=kernel[i];
				kernel[i]=kernel[s-i];
				kernel[s-i]=c;
			}
		}
		
		
	//The filter is a "Direct Form II Transposed"
		//  implementation of the standard difference equation:
		//   a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
		//                         - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
		//  
		// based on Chen Yangquan <elecyq@nus.edu.sg> 1998-11-11
		// http://mechatronics.ece.usu.edu/yqchen/filter.c/
		
		public static double[] filter2(double[]b, double[]a,  double[]x)	{
			double[] y= new double[x.length];
		    int ord=Math.min(a.length, b.length);
		    
		    // initial conditions
		    y[0]=b[0] * x[0];
		    for (int i=1;i<ord;i++){	  
		        for (int j=0;j<=i;j++) {
		            y[i]+=b[j]*x[i-j];
		          //  System.out.print(j+" ");
		        }
		       // System.out.println("\\ "+i);
		        for (int j=1;j<=i;j++) {
		            y[i]-=a[j]*y[i-j];
		        }
		    }
		    /* end of initial part */
		    for (int i=ord;i<x.length;i++) {
		        for (int j=0;j<b.length;j++)
		            y[i]+=b[j]*x[i-j];	        
		        for (int j=1;j<a.length;j++)
		            y[i]-=a[j]*y[i-j];
		    }
		    
		    return y;
		    
		} /* end of filter */
		
		
		public static double[] filter2(double[]b, double a, double[]x)	{
			double[] y= new double[x.length];
		     
		    if (a!=1.0d) {
		    	//System.out.println("normalization");
		    	for (int k=0; k<b.length; k++)
		    		b[k]/=a;
		    }
		   // System.out.println("\nb:");
			//printvector(b);
			
			 // initial conditions
		    y[0]=b[0] * x[0];
		    for (int i=1;i<b.length;i++){	  
		        for (int j=0;j<=i;j++) {
		            y[i]+=b[j]*x[i-j];
		          //  System.out.print(j+" ");
		        }
		      
		    }
		    /* end of initial part */
		    for (int i=b.length;i<x.length;i++) {	    	
		        for (int j=0;j<b.length;j++)  
		            y[i]+=b[j]*x[i-j];	   	          
		    }
		    
		    return y;
		    
		} /* end of filter */
		
		public static double[] filterTF(double[]b, double[]a,  int k){
			k=Math.abs(k);
		    int ord=Math.min(a.length, b.length);
		    double[]x = new double[k*ord];
		    x[0]=1;
		    double[]y= filter2( b,  a,   x);
		    return y;
		}
		
		public static boolean isPow2(int i) {
			if (i%2 == 0) {
				while (i>2) {i=i/2; 
				//System.out.println(i);
				if (i%2==1) return false;
				}
				return true;
			}
			return false;
		}
		
		public static int nfft(int z) {
			return  pow2( nextpow2(z));
		}
		
		public static int nextpow2(int u) {
			int i=0;
			if (u<0) u=-u;
			if (u<=1)
				return 0;
			if ((u&(u-1)) ==0) { i--; }
			while (u>0) {
				u>>=1;
				//System.out.println(u +" "+i);
				i++;
			}
			return i;	
		}
		
		public static int pow2(int n) {
			if (n==0) return 1;
		/*	if ((n & (n - 1)) == 0) {
				return n; // x is already a power-of-two number 
			}*/
			return 2<<(n-1);
		}

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
}
