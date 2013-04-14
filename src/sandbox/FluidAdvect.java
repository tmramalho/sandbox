package sandbox;

public class FluidAdvect implements Runnable {
	private final double[] vxa;
	private final double[] vya;
	private final double[] ain;
	private final double[] aout;
	private final double c;
	private final int rowSize;
	private final int colSize;
	private final boolean[] wall;
	public FluidAdvect(double[] vxa_, double[] vya_, double[] ain_, double[] aout_, double c_, int rowSize_, int colSize_, boolean[] wall_) {
		vxa = vxa_;
		vya = vya_;
		ain = ain_;
		aout = aout_;
		c = c_;
		rowSize = rowSize_;
		colSize = colSize_;
		wall = wall_;
	}
	
	public void run() {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				double vx = vxa[p];
				double vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					double xn = x + vx * c;
					double yn = y + vy * c;
					
					double lb = rowSize - 1.0001;
					double ub = colSize - 1.0001;
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					if(wall[pn]) {
						//if we advect into a wall, we advect instead into ourselves
						xn = x;
						yn = y;
						if(xn < 0) xn = 0;
						else if (xn > lb) xn = lb;
						if(yn < 0) yn = 0;
						else if (yn > ub) yn = ub;
						
						pn = cell((int)xn,(int)yn);
					}
					
					double fx = xn - (int)xn;
					double fy = yn - (int)yn;
					
					double in = ain[p];
					
					double ia = (1.0-fy) * (1.0-fx) * in;
					double ib = (1.0-fy) * fx * in;
					double ic = fy * (1.0-fx) * in;
					double id = fy * fx * in;
					aout[p] -= (ia+ib+ic+id);
					aout[pn] += ia;
					aout[pn+1] += ib;
					aout[pn+rowSize] += ic;
					aout[pn+rowSize+1] += id;
				}
			}
		}
	}
	
	private int cell(int x, int y) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		return x + y * rowSize;
	}
}
