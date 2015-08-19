import java.io.File;

import dsp.*;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ijaux.Constants;
import ijaux.PixLib;
import ijaux.Util;
import ijaux.datatype.Pair;
import ijaux.hypergeom.ComplexCube;
import ijaux.hypergeom.PixelCube;
import ijaux.hypergeom.dsp.FFTD;
import ijaux.hypergeom.index.BaseIndex;
import ijaux.scale.CLineIteratorIP;
import ijaux.scale.IJLineIteratorIP;

/**
 * 
 */

/**
 * @author adminprodanov
 *
 */
public class FFTPlugin_ implements PlugInFilter {

	private String version="1.0";
	final int flags=DOES_8G+
			DOES_8C+
			DOES_16+
			DOES_32+
			NO_CHANGES+
			NO_UNDO;
	private boolean isFloat=false;
	private boolean isRGB=false;
	
	 

	void showAbout() {
		IJ.showMessage("FFT Plugin "+version,
				"The plugin computes FFT"
				);
	}

	@Override
	public int setup(String arg, ImagePlus imp) {
		if (imp==null) return DONE;
		if (arg.equals("about")){
			showAbout();
			return DONE;
		}
		
		if(IJ.versionLessThan("1.48") ) {
			return DONE;
		}
		else {		
			IJ.hideProcessStackDialog=true;
	
			isFloat= (imp.getType()==ImagePlus.GRAY32);
			isRGB = (imp.getType()==ImagePlus.COLOR_RGB);
			
			//makeRedLut();
			return IJ.setupDialog(imp, flags);
		}
	}

	
	
	@Override
	public void run(ImageProcessor ip) {
	 	int width=ip.getWidth();
	 	int height=ip.getHeight();
	 	int k=0;	 	
	 	
	 	if (DSP.isPow2(width) && DSP.isPow2(height)) {
	 		ImageProcessor ipaux;
	 		FloatProcessor fpaux=new FloatProcessor(2*width, height);
			if (!isFloat) 
				ipaux=ip.toFloat(0, null);
			else 
				ipaux=ip;
			
			IJLineIteratorIP<float[]> iter=new IJLineIteratorIP<float[]>(ipaux, 0);
			IJLineIteratorIP<float[]> fiter=new IJLineIteratorIP<float[]>(fpaux, 0);		
			CLineIteratorIP citer=new CLineIteratorIP(fpaux, 1);
			Pair<double[], double[]> ptab=FFTUtil.expTable2(2*width, -1);
			float[] realrow=null, complexrow=null;

			while (iter.hasNext()) {
				realrow=iter.next();
				complexrow=FFTProc.rfftp(realrow, ptab);
				fiter.putLine(complexrow, k);
				k++;
			} // xdir
			k=0;
			while (citer.hasNext()) {
				complexrow=citer.next();
				FFTProc.cfftp(complexrow, true, ptab);
				citer.putLine(complexrow, k);
				k+=2;
			} // ydir
			
			new ImagePlus("FFT", fpaux).show();
	 	} else {
	 		IJ.log("not a power of 2 image");
	 	}
		
	}

	public static float[] complexify (float[] real) {
		int n=2* real.length;
		
		float[] ret = new float[n];
		
		for (int i=0; i<real.length; i++) {
			ret[2*i]=real[i];
		}
		
		return ret;
		
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int width=32,   height=32;
		
		float[] pixels=Util.rampFloat(width*height, 128);
		double[] dpixels=Util.rampDouble(width*height, 128);
		
		FloatProcessor ipaux=new FloatProcessor(width, height, pixels, null);
		
 		FloatProcessor fpaux=new FloatProcessor(2*width, height);
 		//FloatProcessor fpaux2=new FloatProcessor(2*width, height);
		Pair<double[], double[]> ptab=FFTUtil.expTable2(2*width, -1);
		
		IJLineIteratorIP<float[]> iter=new IJLineIteratorIP<float[]>(ipaux, 0);
		IJLineIteratorIP<float[]> fiter=new IJLineIteratorIP<float[]>(fpaux, 0);		
		
		float[] realrow=null, complexrow=null;
		int k=0;
		while (iter.hasNext()) {
			realrow=iter.next();
			complexrow=FFTProc.rfftp(realrow, ptab);
			//TestUtil.printcvector(complexrow);
			//System.out.println("========");
			fiter.putLine(complexrow, k, 0);
			k++;
		} // xdir
		CLineIteratorIP citer=new CLineIteratorIP(fpaux, 1);
		//CLineIteratorIP citer2=new CLineIteratorIP(fpaux2, 1);
		k=0;
		//System.out.println("========");
		//System.out.println("========");
		//System.out.println("========");
		while (citer.hasNext()) {
			//complexrow=citer.getLineFloat(fpaux, k, 1);
			complexrow=citer.next();
			//System.out.println(k+ "  \n>>========");
			//TestUtil.printcvector(complexrow);
			FFTProc.cfftp(complexrow, true, ptab);
			//System.out.println("  \r\n<<========");
			//TestUtil.printcvector(complexrow);
			//citer2.putLine( complexrow, k );
			citer.putLine( complexrow, k );
			k+=2;
			//citer2.fwd();
		} // xdir
		
		PixelCube<Double,BaseIndex> pc= new PixelCube<Double,BaseIndex> (new int[]{width, height}, dpixels);
		pc.setIndexing(Constants.BASE_INDEXING);
		ComplexCube cc=FFTD.fftxd(pc, -1);
		
		PixLib plib=new PixLib();
		
		try {

			File f=new File(args[0]);

			if (f.exists() && f.isDirectory() ) {
				System.setProperty("plugins.dir", args[0]);
				
			} else {
				throw new IllegalArgumentException();
			}
			
			new ImageJ();
			new ImagePlus("orig", ipaux).show();
			
			new ImagePlus("FFT", fpaux).show();
			
			ImagePlus imgp2=plib.imageFrom("FFT 2", cc, Constants.FFT_R, 0);
			imgp2.show();
		}
		catch (Exception ex) {
			IJ.log("plugins.dir misspecified\n");
			ex.printStackTrace();
		}

	}

}
