package test.dsp;
import dsp.FFTProc;
import dsp.FFTProc3;
import dsp.FFTUtil;
import ij.process.FloatProcessor;
import ijaux.Util;
import ijaux.datatype.Pair;
import ijaux.dsp.DSP;
import ijaux.dsp.FFT;
import static dsp.TestUtil.*;

public class RealtestFFT1D
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
	
	static final float[] x6r={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		 
		float[] row4=x.clone();
		float[] row1=x2.clone();
		float[] row2=x3.clone();
		FFTUtil.bit_reverse2(row2);
		System.out.println ("\nBit reverse");
		Util.printFloatArray(row2);
		
		int nfft=FFT.nfft(row1.length);
		FFTUtil.bit_reverse2(row2);
		System.out.println ("\\\\ "+DSP.nextpow2(row1.length)+" nfft "+nfft);
		
 
		System.out.println ("********************************");
		System.out.println ("FFT");
		float[] imag= new float[row1.length];
		float[] imag2= new float[x.length];
		Pair<float[], float[]> carr= 	FFTProc.fftC2C1d(x, imag2, -1, x.length);
		System.out.println ("\nReal part");
		Util.printFloatArray(carr.first);
		System.out.println ("\nImaginary part");
		Util.printFloatArray(carr.second);
		
		System.out.println ("Complex 2 complex interleaved");
		//float[] re=carr.first;
		//float[] im=carr.second;
		FFTProc.cfft(row1, true);
		
		float[] re=new float[ row1.length/2];
		float[] im=new float[ row1.length/2];
		int k=0;
		for (int i=0; i< re.length; i++) {
			re[i]=row1[k];
			im[i]=row1[k+1];
			k+=2;
		}
		
		System.out.println ("\nReal part cfft");
		Util.printFloatArray(re);
		System.out.println ("\nImaginary part cfft");
		Util.printFloatArray(im);
		double r1=corrcoef(re, xr);
		
		double r2=corrcoef(im, xi);
		System.out.println ("cr " +r1 + " test passed: " +(r1==1) );
		System.out.println ("cr " +r2 + " test passed: " +(r2==1) +" combined " +(r1==1  &&  r2==1) );
		
		float[] x3=x2.clone();
		FFTProc.cfftp(x3, true, null);
		k=0;
		for (int i=0; i< re.length; i++) {
			re[i]=x3[k];
			im[i]=x3[k+1];
			k+=2;
		}
		
		System.out.println ("\nReal part cfftp");
		Util.printFloatArray(re);
		System.out.println ("\nImaginary part cfftp");
		Util.printFloatArray(im);
		 r1=corrcoef(re, xr);
		 
		 r2=corrcoef(im, xi);
		System.out.println ("cr " +r1 + " test passed: " +(r1==1) );
		System.out.println ("cr " +r2 + " test passed: " +(r2==1) +" combined " +(r1==1  &&  r2==1) );
	
		/*
		
		carr= FFTProc.fftR2Cp1d(row2, -1, nfft);
		re=carr.first;
		im=carr.second;
		System.out.println ("\nReal part");
		Util.printFloatArray(re);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(im);
		
		 r1=corrcoef(re, xr);
		 
		 r2=corrcoef(im, xi);
		System.out.println ("cr " +r1 + " test passed: " +(r1==1) );
		System.out.println ("cr " +r2 + " test passed: " +(r2==1) +" combined " +(r1==1  &&  r2==1) );
	*/
		FFTProc3.rfft2(row4.clone());
		int n=row4.length;
		Pair<double[], double[]> ptab=FFTUtil.expTable(n, -1);
		//float[] newft=FFTProc.rfftp(row4.clone(), true, ptab);
		float[] newft=FFTProc3.rfftp2(row4.clone(), ptab);
		System.out.println ("\n interleaved");
		Util.printFloatArray(newft);
		Pair<float[], float[]> carr2=complexInline(newft);
		re=carr2.first;
		im=carr2.second;
		
		System.out.println ("\nReal part 3");
		Util.printFloatArray(carr2.first);
		System.out.println ("\nImaginary part 3");
		Util.printFloatArray(carr2.second);
		
		r1=corrcoef(re, xr);
		 
		r2=corrcoef(im, xi);
		System.out.println ("cr " +r1 + " test passed: " +(r1==1) );
		System.out.println ("cr " +r2 + " test passed: " +(r2==1) +" combined " +(r1==1  &&  r2==1) );
	
		
		System.out.println ("********************************");
		System.out.println ("\nIFFT");
		nfft=FFT.nfft(xr.length);
		Pair<float[], float[]> carri= FFTProc.fftC2C1d(xr, xi, +1, nfft);
		
		System.out.println ("\nReal part");
		Util.printFloatArray(carri.first);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(carri.second);
		
		System.out.println ("***********fftC2C1d********************");
		float[] x6i=new float[16]; 
		carri=FFTProc.fftC2C1d(x6r.clone(), x6i, -1, 16);
		
		System.out.println ("\nReal part");
		Util.printFloatArray(carri.first);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(carri.second);
		
		 
		System.out.println ("*************fftC2C4******************");
		
		x6i=new float[16]; 
		FFTProc.fftC2C4(x6r, x6i);
		
		System.out.println ("\nReal part");
		Util.printFloatArray(x6r);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(x6i);
		
		/*x6i=new float[16];
		carri=FFTProc.fftC2C1d(x6r.clone(), x6i, -1, x6r.length);
	 
		System.out.println ("\nReal part");
		Util.printFloatArray(carri.first);
		
		System.out.println ("\nImaginary part");
		Util.printFloatArray(carri.second);*/
		
	}
	 
	
	public static Pair<float[], float[]> complexInline(float[] arr) {
		if (arr.length %2 !=0) throw new IllegalArgumentException ("odd length "+arr.length);
		float[] re=new float[arr.length/2];
		float[] im=new float[arr.length/2];
		int k=0;
		for (int i=0; i<arr.length; i+=2) {
			re[k]=arr[i];
			im[k]=arr[i+1];
			k++;
		}
		return Pair.of(re,im);
		
	}

}
