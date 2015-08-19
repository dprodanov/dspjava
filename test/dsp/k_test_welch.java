package test.dsp;

//import static dsp.TestUtil.corrcoef;
//import ijaux.Util;
//import ijaux.datatype.Pair;
//import ijaux.dsp.DSP;
//import ijaux.dsp.FFT;
import static dsp.TestUtil.corrcoef;

import java.util.Arrays;

import ijaux.Util;
import ijaux.datatype.Pair;
import ijaux.dsp.DSP;
import ijaux.dsp.FFT;
import dsp.FFTProc;
//import dsp.kFFTProc;

import static java.lang.Math.*;

public class k_test_welch
{

	//static final float[] x1 ={1, 2,	3,	9,	8,	5,	1,	2};
	static final double[] x1 ={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	static final float[] x2 ={0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0};
	static final float[] g1  ={1,2,3,4,5,6,7,8};
	static final float[] XX ={9,12,-4.7320f,-1.2679f,-1.2679f,-4.7320f,9,12};
	
	
	public static double[][] Segmentation (double[] g, int L, double D){
	    /*
	    x=[1,2,3,4,5,6,7,8,9,10,11,12]
	    x=range(100)
	    L=4
	    D=0.5 -> 0.25, 0.125, 0.0625 
	    overlap = 'regular' or 'circular'
	    */
		if ((L%2!=0) && (L!=1))      
			throw new IllegalArgumentException("L must be power of two" + L);
		
		
		if (((1/D)%2!=0) && ((1/D)!=1))    
		    throw new IllegalArgumentException("D must be power of two " + D);
		
		
		int offset = 0;
		// M - length of one segment
		int M = g.length/L;
		int mD=(int)(D*M);
		//System.out.println ("mD: "+mD);
	    if (mD==0)
	    	mD=1;
	    //	throw new IllegalArgumentException("mD is equal zero");
		// K - number of segments
		//int K = (int) (1/D) * (L-1) + 1;
	    int K = (int) ((L-1)*(M/mD))+7;
	    
	    
		
	    K=0;
	    while ((M+offset)<=g.length) {
	    	offset= offset+mD;
	    	K=K+1;
	   }
	    
	    
	    double[][] S = new double [K][M];
	    
	    offset = 0;
		for (int i=0; i<K; i++){
			//System.out.println ("mD: "+(M+offset));
			double[] el = Arrays.copyOfRange(g,offset, M+offset);
			S[i]=el;
			offset= offset+mD;
			
			//if ((M+offset)>x.length){
			//	System.out.println ("M+off: "+(M+offset)+"off: " );
			//	break;
			//}
		}
		return S;
	}

	
	public static void main(String[] args) {
		
		// Choose the input array here:
		double[] g = x1.clone();
		//int L   = 2;
		//float D = 0.5f;
		int L = 2;
		int D = 4;
		/***********************************************/
		double[][] S  = Segmentation(g,L,D);
		/*
		 float[][] S  = kFFTProc.Segmentation(g,L,D);
		float[]   W  = kFFTProc.HanningWindow (S[0].length);
		float[][] F  = kFFTProc.OperateOnSegments (S,W);
		float[]  SM  = kFFTProc.SquaredMagnitude (F[0]);
		float   SCF  = kFFTProc.ComputeScaleFactor (W, g.length);
		float[] SCSM = kFFTProc.ApplyScale (SM, SCF);
		float[]   R  = kFFTProc.Apply10LOG10 (SCSM);
		float[] FinalResult = kFFTProc.kWelch ( g, L, D);
		float[] ZP =   kFFTProc.ZeroPadding (g, 0);
		float[] FRZP = kFFTProc.kWelch ( ZP, L, D);
		*/
		/****  PRINT THE RESULTS *****/
		//System.out.println ("---Scale factor: "+S[14]);
		for (int i=0; i<S.length; i++){
			Util.printDoubleArray(S[i]);
		}
		/*
		System.out.println ("---Hanning function: "+W.length);
		Util.printFloatArray(W);
		System.out.println ("---FFT array: ");
		for (int i=0; i<3; i++){
			Util.printFloatArray(F[i]);
		}
		System.out.println ("---SM: ");
		Util.printFloatArray(SM);
		System.out.println ("---Scale factor: ");
		System.out.println (""+ SCF);
		System.out.println ("---Scaled array: ");
		Util.printFloatArray(SCSM);
		System.out.println ("---one result: ");
		Util.printFloatArray(R);
		System.out.println ("---final result: ");
		Util.printFloatArray(FinalResult);
		//System.out.println ("---Simple test: "+log10(exp(1)));
		//System.out.println ("---Zero Padding: ");
		//Util.printFloatArray(ZP);
		System.out.println ("---Final result after Zero Padding: ");
		Util.printFloatArray(FRZP);
		*/
		/****  END OF PRINT *****/
	}
	
}