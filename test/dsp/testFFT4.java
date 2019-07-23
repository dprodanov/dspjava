package test.dsp;
import dsp.FFTProc;
import dsp.FFTUtil;
import ijaux.TestUtil;
import ijaux.Util;
import ijaux.datatype.Pair;

public class testFFT4
 {
	/*
	 *  FFT of [ 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 ]
	 *  in Matlab
	 */
	
	/* re
	 * [120,-8,-8.0,-8.000000000000002,-7.999999999999999,-7.999999999999998,-
7.999999999999997,-7.999999999999996,-8,-8,-8.0,-7.999999999999998,-8.0,-
8.000000000000002,-8.000000000000004,-8.000000000000004]

[120.0,-8.0,-8.0,-8.000000000000002,-7.999999999999999,-7.999999999999998,-
7.999999999999997,-7.999999999999996,-8.0,-8.0,-8.0,-7.999999999999998,-8.0,-
8.000000000000002,-8.000000000000004,-8.000000000000004]
	 */
	static float[] xr={
			120.0f,	-8.0f,	-8.0f,	-8.0f,
			-8.0f,	-8.0f,	-8.0f,	-8.0f,
			-8.0f,	-8.0f,	-8.0f,	-8.0f,
			-8.0f,	-8.0f,	-8.0f,	-8.0f
			 

	};
	static float[] xi={
			0.0f,				40.2187159370068f,
			19.3137084989848f,	11.9728461013239f,
			8.0f,				5.34542910335439f,
			3.31370849898476f,	1.59129893903727f,
			0.0f,				-1.59129893903727f,
			-3.31370849898476f,	-5.34542910335439f,
			-8.0f,				-11.9728461013239f,
			-19.3137084989848f,	-40.2187159370068f

	};

	/* im
	 * [0, 40.21871593700678,19.31370849898476,11.97284610132391,8.0,
5.345429103354393,3.313708498984761,1.591298939037266,0,-1.591298939037262,
-3.313708498984761,-5.345429103354393,-8.0,-11.97284610132391,-19.31370849898476,
-40.21871593700678]

[0.0,40.21871593700678,19.31370849898476,11.97284610132391,8.0,
5.345429103354393,3.313708498984761,1.591298939037266,0.0,-1.591298939037262,-
3.313708498984761,-5.345429103354393,-8.0,-11.97284610132391,-19.31370849898476,-
40.21871593700678]
	 */
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		float [] arr=Util.rampFloat(16, 16);
		Util.printFloatArray(arr); 
		float[] cplx=complexify(  arr);
		
		System.out.println ("**************complexification******************");
		System.out.println ("interleaved");
		Util.printFloatArray(cplx); 

		
		System.out.println ("********************************");
	 
		Pair<double[], double[]> ptab=FFTUtil.expTable2(cplx.length, -1);
		System.out.println ("CFFT");
		 FFTProc.cfftp(cplx, true, ptab);
		
		 Util.printFloatArray(cplx); 
			System.out.println ("re");
		float[] re2=getRe(cplx); 
		Util.printFloatArray(re2); 
		Util.printFloatArray(xr); 
		System.out.println("corr re: " +TestUtil.corrcoef(re2, xr));

		float[] im2=getIm(cplx); 
		System.out.println ("im");
		Util.printFloatArray(im2); 
		Util.printFloatArray(xi); 
		System.out.println("corr im: " +TestUtil.corrcoef(im2, xi));
		
		System.out.println ("ICFFT");
		ptab=FFTUtil.expTable2(cplx.length, +1);
		 FFTProc.cfftp(cplx, false,ptab);
		
		 
		 float[] re3=getRe(cplx); 
		 System.out.println ("re");
		 Util.printFloatArray(re3); 
		
		float[] im3=getIm(cplx); 
		System.out.println ("im");
		Util.printFloatArray(im3); 
		
		ptab=FFTUtil.expTable(16, +1);
		System.out.println ("expTable");
		Util.printDoubleArray(ptab.first); 
		Util.printDoubleArray(ptab.second); 
		Pair<double[], double[]> ptab2=FFTUtil.expTable2(16, +1);
		System.out.println ("expTable2");
		Util.printDoubleArray(ptab2.first); 
		Util.printDoubleArray(ptab2.second); 
	}
	
	static float[] complexify(float[] arr) {
		float[] ret=new float[ 2* arr.length];
		
		for (int i=0; i<arr.length; i++) 
			ret [2*i]=arr[i];
		
		return ret;
		
	}
	
	static float[] getRe(float[] arr) {
		float[] ret=new float[ arr.length/2];
		
		for (int i=0; i<arr.length; i+=2) 
			ret [i/2]=arr[i];
		
		return ret;
		
	}
	
	static float[] getIm(float[] arr) {
		float[] ret=new float[ arr.length/2];
		
		for (int i=1; i<arr.length; i+=2) 
			ret [i/2]=arr[i];
		
		return ret;
		
	}
	
}
