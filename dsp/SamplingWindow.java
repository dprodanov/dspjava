package dsp;

import ijaux.datatype.UnsupportedTypeException;
import ijaux.datatype.oper.DoubleOp;
import ijaux.datatype.oper.Op;
import static java.lang.Math.*;

/*
 * Sampling windows
 * 
 */
public abstract class SamplingWindow {
	
	protected double [] wnd;
	
	DoubleOp op;
	
	protected int n=0;
	
	public static final double TWOPI = 2*Math.PI;
	
	public static SamplingWindow createWindow(WindowTypes type, int n) {		
		switch (type ) {
			case HANNING: return new HanningWindow(n);
			case HAMMING: return new HammingWindow(n);
			case LANCZOS: return new LanczosWindow(n);
			case GAUSSIAN: return new GaussianWindow(n, 1.0);
		}
		
		return null;
			
	}
	
	public SamplingWindow(int n) {
		if (n==1) 
			throw new IllegalArgumentException("n>1 needed");
		
		this.n=n;
		wnd=new double[n];
		 try {
			op= (DoubleOp) Op.get(double.class);
		} catch (UnsupportedTypeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		compute();
	}
	
	

	public double[] getWindow() {
		return wnd;
	}
	
	public void apply(double[] array) {
		op.mult(array, wnd);
	}
	
	public abstract void compute ();
	
	public static double sinc(double x) {
		if (x==0) return 1.0d;
		final double px=PI*x;
		return sin(px)/px;
	}
	
	@Override
	public String toString() {
		StringBuffer sb=new StringBuffer();
		sb.append("[ ");
		for (double d: wnd) {
			sb.append(d+" ");
		}
		sb.append("]");
		return sb.toString();
	}
	
}

/*
 * Hanning window
 */
class HanningWindow extends SamplingWindow {

	public HanningWindow(int n) {
		super(n);
	}

	@Override
	public void compute() {
		for (int i=0; i<n; i++){
			wnd[i]=0.5*(1-cos(TWOPI* i/ (n-1)));
		}
		
	}

} //


/*
 * Hamming window
 */
class HammingWindow extends SamplingWindow {

	public HammingWindow(int n) {
		super(n);
	}
	
	// https://ccrma.stanford.edu/~jos/sasp/Hamming_Window.html
	final double alpha=25.0/46.0, beta=21.0/92.0;
	
	@Override
	public void compute() {
		for (int i=0; i<n; i++){
			wnd[i]= alpha - beta*cos(TWOPI* i/ (n-1));
		}
		
	}
	


} //

/*
 * Lanczos window
 */
class LanczosWindow extends SamplingWindow {

	public LanczosWindow(int n) {
		super(n);
	}

	@Override
	public void compute() {
		for (int i=1; i<n-1; i++){
			double arg= 2* i;
			//System.out.println(arg/ (n-1) -1);
			arg=arg/ (n-1) -1;
			wnd[i]= sinc(arg);			
		}
		
	}

} //

/*
 * Gaussian window
 */
class GaussianWindow extends SamplingWindow {
	private boolean computed =false;
	public GaussianWindow(int n, double s) {
		super(n);
		this.s=s;
		//System.out.println("gaussian "+s);
		computed=true;
		compute();
	}
	
	double s=0.5;

	@Override
	public void compute() {
		if (! computed)
			return;
		
		//System.out.println("gaussian 1 "+s);
		for (int i=0; i<n; i++){
			final double x=(2.0* (double)i/(double) (n-1) -1)/s;
			//System.out.println("gaussian "+x);
			wnd[i]=exp( -x*x);
		}
		
	}

} //



