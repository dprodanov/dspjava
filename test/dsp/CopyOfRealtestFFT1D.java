package test.dsp;
import dsp.DSP;
import dsp.FFTUtil;
import ij.process.FloatProcessor;
 
import ijaux.datatype.Pair;
import static dsp.TestUtil.*;

public class CopyOfRealtestFFT1D
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
	
	// FFT by rows
	static float[] xr_re={15,    -2,    -7,    -2,
		    16,     7,     2,     7,
		    26,     2,     6,     2};
	
	static float[] xr_im={0,     7,     0,    -7,
		     0,    -3,     0,     3,
		     0,    -6,     0,     6};
	

	
/*	static float[] arr={1,  2,  3,  9,
			8,  5,  1, 2,
			9 , 8 , 7,  2};*/
	
	static final float[] x={1,2,3,9,8,5,1,2};
	static final float[] x2={1,0, 2,0, 3,0, 9,0, 8,0, 5,0, 1,0, 2,0};
	static final float[] x3={1,2,3,4,5,6,7,8};
	static final float[] xr={31,	-14.0710678118655f,	5,	0.0710678118654755f,	
			-5,	0.0710678118654755f,	5,	-14.0710678118655f	};
	
	static final float[] xi={0,	-4.82842712474619f,	4,	-0.828427124746190f,	
		0,	0.828427124746190f,	-4,	4.82842712474619f};
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		System.out.println ("is 2^x " +4+" " +DSP.isPow2(4));
		System.out.println ("is 2^x " +6+" " +DSP.isPow2(6));
		System.out.println ("is 2^x " +17+" " +DSP.isPow2(17));
		System.out.println ("is 2^x " +32+" " +DSP.isPow2(32));
	
		System.out.println ("Twiddele table");
		System.out.println ("********************************");
		
		System.out.println ("FFTProc.expTable2");
		int n=16;
		Pair<double[], double[]> ptab=FFTUtil.expTable2(n, -1);
		double[] cosTable =ptab.first;
 		double[] sinTable =ptab.second;
		int k=cosTable.length;
		for (int j=0; j<k; j++) {
			System.out.println(j +" cos "+ cosTable[j]+" sin " + sinTable[j]);
		}
		System.out.println ("********************************");
		System.out.println ("FFTProc.expTable");
		ptab=FFTUtil.expTable(n, -1);
		cosTable =ptab.first;
 		sinTable =ptab.second;
	
		for (int j=0; j<k; j++) {
			System.out.println(j +" cos "+ cosTable[j]+" sin " + sinTable[j]);
		}
		///////////////////////
		System.out.println ("********************************");
		System.out.println ("FFTProc.iexpTable");
		double[] wt=FFTUtil.iexpTable(n, -1);
		
		for (int j=0; j<wt.length; j+=2) {
			System.out.println( (j/2) +" cos "+ wt[j]+" sin " + wt[j+1]);
		}
		
		System.out.println ("FFTProc.iexpTable2");
		wt=FFTUtil.iexpTable2(n, -1);
		
		for (int j=0; j<wt.length; j+=2) {
			System.out.println( (j/2) +" cos "+ wt[j]+" sin " + wt[j+1]);
		}
		
	}
	 
	
	
}
