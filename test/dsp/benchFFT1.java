package test.dsp;
import dsp.FFTProc;
import dsp.FFTProc3;
import dsp.FFTUtil;
import ijaux.Util;
import ijaux.datatype.Pair;
import ijaux.dsp.FFT;

public class benchFFT1 {

	 private static int[] sizes1D = new int[] {32768, 65536, 131072, 262144, 524288,1048576, 2097152};
	 private static int[] sizes1D4 = new int[] {1024, 4096, 16384, 65536, 262144, 1048576};

	 private static int niter = 2; // 200
	/**
	 * @param args
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws InterruptedException {
		System.gc();
		System.gc();
	
System.out.println("******************");
		
		System.out.println("FFTProc.cfft>>");
		
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k];
			 
			float[] x=(float[]) Util.rand(n, float.class);
			//float[] row1=new float[x.length*2];
			 
			//int nfft=FFT.nfft(row1.length);
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				FFTProc.cfft(x,   true);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
		
		System.out.println("******************");
		
		
		System.out.println("FFTProc.rfft2>>");
		//Pair<double[],double[]> ptab=null;
		
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k]/2;
			float[] x=(float[]) Util.rand(n, float.class);
 
	 		//int nfft=FFT.nfft(row1.length);
			
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				FFTProc3.rfft2(x);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
		
		System.out.println("******************");
		
		System.out.println("FFTProc.fftR2C1d>>");
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k]/2;
			float[] row1=(float[]) Util.rand(n, float.class);
			int nfft=FFT.nfft(row1.length);
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				
				Pair<float[], float[]> 	carr= FFTProc.fftR2C1d(row1, -1, nfft);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
		
		System.out.println("FFTProc.fftC2C4>>");
		for (int k=0; k<sizes1D4.length; k++) {
			int n=sizes1D4[k];
			float[] row1=(float[]) Util.rand(n, float.class);
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				float[] im=new float[row1.length];
				FFTProc.fftC2C4(row1, im  );
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
		
		System.out.println("******************");
		System.out.println("FFTProc.fftR2Cp1d>>");
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k];
			float[] row1=(float[]) Util.rand(n, float.class);
			int nfft=FFT.nfft(row1.length);
			long time=-System.nanoTime();		
			for (int i=0; i<niter; i++) {
				FFTProc.fftR2Cp1d(row1, -1, nfft);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
	
		System.out.println("******************");
		/*
		System.out.println("FFTProc.cfftp>>");
		
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k];
			Pair<double[],double[]> ptab=FFTProc.expTable(2*n, -1);
			float[] x=(float[]) Util.rand(n, float.class);
			float[] row1=new float[x.length*2];
			 
			//int nfft=FFT.nfft(row1.length);
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				FFTProc.cfftp(row1,   true, ptab);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
	*/
		System.out.println("******************");
		
		
		System.out.println("FFTProc.rfftp2>>");
		//Pair<double[],double[]> ptab=null;
		
		for (int k=0; k<sizes1D.length; k++) {
			int n=sizes1D[k];
			float[] x=(float[]) Util.rand(n, float.class);
			Pair<double[],double[]> ptab=FFTUtil.expTable(2*n, -1);
	 		//int nfft=FFT.nfft(row1.length);
			
			long time=-System.nanoTime();	
			for (int i=0; i<niter; i++) {
				FFTProc3.rfftp2(x,  ptab);
			}
			time+=System.nanoTime();
			System.out.println("size: "+n+"\t execution time: " + String.format("%.4f", time / 1000.0/niter ) + " usec");
			
		}
		
		System.out.println("******************");
		
	}

}
