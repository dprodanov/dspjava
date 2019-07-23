package test.dsp;
import dsp.FFTProc;
import ijaux.Util;
import ijaux.datatype.Pair;
import dsp.DSP;


public class testFFT2 {
	// FFT by colls
	
	static float[] xc_re={120,	-8,	-8,	-8,	
						   -8,	-8,	-8,	-8,	
						   -8,  -8,	-8,	-8,	
						   -8,	-8,	-8,	-8};

	static float[] xc_im={0,	40.2187159370068f,	19.3137084989848f,	11.9728461013239f,	
						  8,	5.34542910335439f,	3.31370849898476f,	1.59129893903727f,
						  0,	-1.59129893903727f,	-3.31370849898476f,	-5.34542910335439f,
						 -8,	-11.9728461013239f,	-19.3137084989848f,	-40.2187159370068f};
	 
	static float[] arr=Util.rampFloat(16, 16);
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Util.printFloatArray(arr);
		
		int nfft=DSP.nfft(arr.length);
		System.out.println ("nfft :" +nfft);
		final Pair<float[], float[]> carr= FFTProc.fftR2C1d(arr,  -1, nfft);
		float[] re=carr.first;
		float[] im=carr.second;
	/*	System.out.println ("\nReal part");
		Util.printFloatArray(re);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(im);
		*/
		printArr2(re, im);
		FFTProc pr =new FFTProc(nfft);
		int N=nfft;
		float[] xarr=new float[2*N];
		System.arraycopy(arr, 0, xarr, 0, arr.length);	
				
		//pr.rfft(xarr, N, true);
		System.out.println ("\n nested:");
		Util.printFloatArray(xarr);
		
		/*final Pair<float[], float[]> carr2= FFTProc.fftR2Cp1d(arr,  -1, nfft);
		float[] re2=carr2.first;
		float[] im2=carr2.second;
		System.out.println ("\nReal part");
		Util.printFloatArray(re2);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(im2);*/
	 	 
	}
	
	static void printArr2(float[] a, float[] b) {
		System.out.print("\n[\n");
		for (int i=0; i< a.length; i++) {
			System.out.print("("+a[i]+", "+ b[i] +")\n");
		}
		System.out.print("\n]\n");
	}

}
