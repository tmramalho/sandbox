package sandbox;
import processing.core.PApplet;

public class Fluid {
	PApplet parent;
	private int cellw;
	private int cellh;
	private int rowSize;
	private int colSize;
	private int tsize;
	private double[] dens0;
	private double[] dens1;
	private double[] vx0;
	private double[] vx1;
	private double[] vx2;
	private double[] vy0;
	private double[] vy1;
	private double[] vy2;
	private double[] req;
	private double[] temp0;
	private double[] temp1;
	private double[] indReq;
	private int[] sourceCell;
	private boolean[] wall;
	private int[] numParticles;
	private double baseDens;
	private double addVel;
	
	public double[] getVx0() {
		return vx0;
	}

	public double[] getVy0() {
		return vy0;
	}
	
	public int[] getNumParticles() {
		return numParticles;
	}
	
	public Fluid(int cw, int ch, int w, int h, PApplet p, double bd) {
		rowSize = cw;
		colSize = ch;
		parent = p;
		cellw = Math.round(w / (rowSize - 2));
		cellh = Math.round(h / (colSize - 2));
		
		tsize = rowSize*colSize;
		dens0 = new double[tsize];
		dens1 = new double[tsize];
		vx0 = new double[tsize];
		vx1 = new double[tsize];
		vx2 = new double[tsize];
		vy0 = new double[tsize];
		vy1 = new double[tsize];
		vy2 = new double[tsize];
		req = new double[tsize];
		temp0 = new double[tsize];
		temp1 = new double[tsize];
		indReq = new double[4*tsize];
		sourceCell = new int[tsize];
		wall = new boolean[tsize];
		numParticles = new int[tsize];
		baseDens = bd;
		addVel = 10;
		
		for(int i = 0; i<tsize; i++) {
			dens0[i] = dens1[i] = bd;
		}
	}
	
	public void update(double a) {
		clearWall();
		udiffuse(dens0, dens1, 0.1*a);
		udiffuse(temp0, temp1, 0.1*a);
		diffuse(vx0, vx2, 0.01*a);
		diffuse(vy0, vy2, 0.01*a);
		
		pressure(dens1, vx2, vy2, 0.05*a);
		temp(temp1, vx2, vy2, 0.1*a);
		
		edgeVelocities(vx2, vy2);
		
		vconf(vx2, vy2, vx0, vy0, 0.01*a);
		
		friction(vx2, vy2, vx1, vy1, 0.1*a, 0.1*a, 0.01*a);
		
		copy(vx1,vx0);
		copy(vy1,vy0);
		
		Runnable fax = new FluidAdvect(vx0, vy0, vx1, vx2, a, colSize, rowSize, wall);
		Runnable fay = new FluidAdvect(vx0, vy0, vy1, vy2, a, colSize, rowSize, wall);
		
		Thread t1 = new Thread(fax);
		Thread t2 = new Thread(fay);
		
		t1.start();
		t2.start();
		
		try{
			t1.join();
		} catch(InterruptedException e) {
			t1.interrupt();
		}
		try{
			t2.join();
		} catch(InterruptedException e) {
			t2.interrupt();
		}
		
		/*advect(vx0, vy0, vx1, vx2, a);
		advect(vx0, vy0, vy1, vy2, a);*/
		
		edgeVelocities(vx2, vy2);
		
		copy(vx2, vx1);
		copy(vy2, vy1);
		
		Runnable fbx = new FluidBadvect(vx0, vy0, vx1, vx2, a, colSize, rowSize, wall);
		Runnable fby = new FluidBadvect(vx0, vy0, vy1, vy2, a, colSize, rowSize, wall);
		
		Thread t3 = new Thread(fbx);
		Thread t4 = new Thread(fby);
		
		t3.start();
		t4.start();
		
		try{
			t3.join();
		} catch(InterruptedException e) {
			t3.interrupt();
		}
		try{
			t4.join();
		} catch(InterruptedException e) {
			t4.interrupt();
		}
		
		/*badvect(vx0, vy0, vx1, vx2, a);
		badvect(vx0, vy0, vy1, vy2, a);*/
		
		edgeVelocities(vx2, vy2);
		
		copy(vx2, vx0);
		copy(vy2, vy0);
		
		copy(dens1, dens0);
		advect(vx0, vy0, dens1, dens0, a*2);
		copy(dens0, dens1);
		buadvect(vx0, vy0, dens1, dens0, a*2);
		copy(dens0, dens1);
		
		copy(temp1, temp0);
		advect(vx0, vy0, temp1, temp0, a*2);
		copy(temp0, temp1);
		buadvect(vx0, vy0, temp1, temp0, a*2);
		copy(temp0, temp1);
	}
	
	public void draw() {
		//leave 1 cell border around the screen
		for(int x = 1; x < rowSize - 1; x++) {
			for(int y = 1; y < colSize - 1; y++) {
				int p = cell(x,y);
				int dens = (int)dens0[p]*2;
				int wallc = 0;
				if(wall[p]) wallc = 255;
				parent.fill(dens, dens + wallc, dens);
				parent.rect((x-1)*cellw, (y-1)*cellh, cellw, cellh);
			}
		}
	}
	
	public void addDens(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		int dens = (int) dens0[p];
		dens += 10*baseDens;
		dens0[p] = (double) dens;
	}
	
	public void addDens(int mouseX, int mouseY, double bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		double dens = dens0[p];
		dens += 10*bd;
		dens0[p] = dens;
	}
	
	public void addTemp(int mouseX, int mouseY, double bd) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		double temp = temp0[p];
		temp += bd;
		temp0[p] = temp;
	}
	
	public void addVel(int mouseX, int mouseY, double angle) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		double nextX = addVel*Math.cos(angle);
		double nextY = addVel*Math.sin(angle);
		vx0[p] = nextX;
		vy0[p] = nextY;
	}
	
	public void updateWalls() {
		int th = (int) Math.floor(0.5 * cellw * cellh);
		for (int p = 0; p < tsize; p++) {
			if(numParticles[p] > th && !wall[p]) wall[p] = true;
			else if(numParticles[p] < th && wall[p]) wall[p] = false;
			numParticles[p] = 0;
		}
	}
	
	public void addWall(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		wall[p] = true;
	}
	
	public void remWall(int mouseX, int mouseY) {
		int posx = (int)(mouseX / cellw) + 1;
		int posy = (int)(mouseY / cellh) + 1;
		int p = cell(posx, posy);
		if(wall[p]) {
			wall[p] = false;
			double acc = 0;
			double div = 0;
			int pxu = ncell(posx, posy, posx+1,posy);
			if(pxu != p) {
				acc += dens0[pxu];
				div += 1.0;
			}
			int pxd = ncell(posx, posy, posx-1,posy);
			if(pxd != p) {
				acc += dens0[pxd];
				div += 1.0;
			}
			int pyu = ncell(posx, posy, posx,posy+1);
			if(pyu != p) {
				acc += dens0[pyu];
				div += 1.0;
			}
			int pyd = ncell(posx, posy, posx,posy-1);
			if(pyd != p) {
				acc += dens0[pyd];
				div += 1.0;
			}
			if(div > 0) dens0[p] = acc / div;
			else dens0[p] = 10;
		}
	}
	
	public void copy(double[] ain, double[] aout) {
		for(int i = 0; i < tsize; i++) {
			aout[i] = ain[i];
		}
	}
	
	public void clearWall() {
		for(int i = 0; i < tsize; i++) {
			if(wall[i]) {
				dens0[i] = 0;
				dens1[i] = 0;
				vx0[i] = 0;
				vx1[i] = 0;
				vx2[i] = 0;
				vy0[i] = 0;
				vy1[i] = 0;
				vy2[i] = 0;
			}
		}
	}
	
	public void reset() {
		for(int i = 0; i < tsize; i++) {
			wall[i] = false;
			dens0[i] = baseDens;
			dens1[i] = baseDens;
			temp0[i] = 0;
			temp1[i] = 0;
			vx0[i] = 0;
			vx1[i] = 0;
			vx2[i] = 0;
			vy0[i] = 0;
			vy1[i] = 0;
			vy2[i] = 0;
		}
	}
	
	public void friction(double[] vxa, double[] vya, double[] vxb, double[] vyb, double a, double b, double c) {
		for(int i = 0; i < tsize; i++) {
			double len = Math.sqrt(vxa[i]*vxa[i]+vxb[i]*vxb[i]);
			if(len == 0) continue;
			double flen = len - len*len*a - len*b - c;
			if(flen < 0) flen = 0;
			double nlen = flen/len;
			vxb[i] = vxa[i] * nlen;
			vyb[i] = vya[i] * nlen;
		}
	}
	
	public void temp(double[] dens, double[] vx, double[] vy, double a) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				int pyd = cell(x,y-1);
				double ay = Math.abs(dens[p]) - Math.abs(dens[pyd]);
				vy[p] -= a*ay;
			}
		}
	}
	
	public void fall(double[] dens, double[] vx, double[] vy, double a) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				double ay = Math.abs(dens[p]);
				vy[p] += a*ay;
			}
		}
	}
	
	public void edgeVelocities(double[] vx, double[] vy) {
		for(int x = 0; x < rowSize; x++) {
			int i = cell(x,0);
			if(vy[i] < 0) vy[i] = -vy[i];
			int j = cell(x,colSize-1);
			if(vy[j] > 0) vy[j] = -vy[j];
		}
		for(int y = 0; y < colSize; y++) {
			int i = cell(0,y);
			if(vx[i] < 0) vx[i] = -vx[i];
			int j = cell(rowSize-1,y);
			if(vx[j] > 0) vx[j] = -vx[j];
		}
		
		/***********************
		 * velocities for walls
		 ***********************/
		
		/*for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) {
					int pxu = ncell(x, y, x+1,y);
					int pxd = ncell(x, y, x-1,y);
					int pyu = ncell(x, y, x,y+1);
					int pyd = ncell(x, y, x,y-1);
					if(vx[pxu] < 0 && !wall[pxu]) vx[pxu] = -vx[pxu];
					if(vx[pxd] > 0 && !wall[pxd]) vx[pxd] = -vx[pxd];
					if(vy[pyu] < 0 && !wall[pyu]) vy[pyu] = -vy[pyu];
					if(vy[pyd] > 0 && !wall[pyd]) vy[pyd] = -vy[pyd];
				}
			}
		}*/
	}
	
	public void pressure(double[] dens, double[] vx, double[] vy, double a) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				double ax = Math.abs(dens[pxd]) - Math.abs(dens[pxu]);
				double ay = Math.abs(dens[pyd]) - Math.abs(dens[pyu]);
				vx[p] += a*ax;
				vy[p] += a*ay;
			}
		}
	}
	
	public void vconf(double[] vxa, double[] vya, double[] w, double[] wa, double c) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				double cx = 0.5*(vxa[pyu] - vxa[pyd]);
				double cy = 0.5*(vya[pxu] - vya[pxd]);
				w[p] = cy - cx;
				wa[p] = Math.abs(cy - cx);
			}
		}
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				double dwx = 0.5*(wa[pxu] - wa[pxd]);
				double dwy = 0.5*(wa[pyu] - wa[pyd]);
				double dn = 1;
				double n = Math.sqrt(dwx*dwx+dwy*dwy);
				if(n != 0) {
					dn = 1.0/n;
				}
				vxa[p] += c*(dwy * dn * w[p]);
				vya[p] -= c*(dwx * dn * w[p]);
			}
		}
		
	}
	
	public void udiffuse(double[] ain, double[] aout, double k) {
		double ki = 1/k;
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int i = cell(x,y);
				if(wall[i]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				double d = k*(ain[pxu] + ain[pxd] + ain[pyu] + ain[pyd] - (4.0 - ki) * ain[i]);
				//aout[i] = d;
				aout[i] = d < 0 ? 0 : d;
			}
		}
	}
	
	public void diffuse(double[] ain, double[] aout, double k) {
		double ki = 1/k;
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int i = cell(x,y);
				if(wall[i]) continue;
				int pxu = ncell(x, y, x+1,y);
				int pxd = ncell(x, y, x-1,y);
				int pyu = ncell(x, y, x,y+1);
				int pyd = ncell(x, y, x,y-1);
				double d = k*(ain[pxu] + ain[pxd] + ain[pyu] + ain[pyd] - (4.0 - ki) * ain[i]);
				aout[i] = d;
				//aout[i] = d < 0 ? 0 : d;
			}
		}
	}
	
	
	public void advect(double[] vxa, double[] vya, double[] ain, double[] aout, double c) {
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
	
	public void badvect(double[] vxa, double[] vya, double[] ain, double[] aout, double c) {
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				double vx = vxa[p];
				double vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					double xn = x - vx * c;
					double yn = y - vy * c;
					
					double lb = rowSize - 1.0001;
					double ub = colSize - 1.0001;
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					double fx = xn - (int)xn;
					double fy = yn - (int)yn;
					
					double ia = (1.0-fy) * (1.0-fx) * ain[pn];
					double ib = (1.0-fy) * fx * ain[pn+1];
					double ic = fy * (1.0-fx) * ain[pn+rowSize];
					double id = fy * fx * ain[pn+rowSize+1];
					
					aout[p] += (ia+ib+ic+id);
					aout[pn] -= ia;
					aout[pn+1] -= ib;
					aout[pn+rowSize] -= ic;
					aout[pn+rowSize+1] -= id;
				}
			}
		}
	}
	
	public void buadvect(double[] vxa, double[] vya, double[] ain, double[] aout, double c) {
		for(int p = 0; p < tsize; p++) req[p]= 0;
		
		for(int x = 0; x < rowSize; x++) {
			for(int y = 0; y < colSize; y++) {
				int p = cell(x,y);
				if(wall[p]) continue;
				double vx = vxa[p];
				double vy = vya[p];
				
				if(vx != 0.0 || vy != 0.0) {
					
					double xn = x - vx * c;
					double yn = y - vy * c;
					
					double lb = rowSize - 1.0001;
					double ub = colSize - 1.0001;
					
					if(xn < 0) xn = 0;
					else if (xn > lb) xn = lb;
					if(yn < 0) yn = 0;
					else if (yn > ub) yn = ub;
					
					int pn = cell((int)xn,(int)yn);
					
					double fx = xn - (int)xn;
					double fy = yn - (int)yn;
					
					double ia = (1.0-fy) * (1.0-fx);
					double ib = (1.0-fy) * fx;
					double ic = fy * (1.0-fx);
					double id = fy * fx;
					
					//save sources and what each loses
					sourceCell[p] = pn;
					
					indReq[4*p] = ia;
					indReq[4*p+1] = ib;
					indReq[4*p+2] = ic;
					indReq[4*p+3] = id;
					
					//accumulate how much each cell loses in total
					req[pn] += ia;
					req[pn+1] += ib;
					req[pn+rowSize] += ic;
					req[pn+rowSize+1] += id;
				} else {
					sourceCell[p] = -1;
				}
			}
		}
		
		for(int p = 0; p < tsize; p++) {
			if(wall[p]) continue;
			int pn = sourceCell[p];
			if(pn != -1) {
				//recover previous data
				double ia = indReq[4*p];
				double ib = indReq[4*p+1];
				double ic = indReq[4*p+2];
				double id = indReq[4*p+3];
				
				//get total fractions and rescale requests
				double fa = req[pn];
				double fb = req[pn+1];
				double fc = req[pn+rowSize];
				double fd = req[pn+rowSize+1];
				
				if (fa<1.0f) fa = 1.0f;
				if (fb<1.0f) fb = 1.0f;
				if (fc<1.0f) fc = 1.0f;
				if (fd<1.0f) fd = 1.0f;
				
				ia = ia * ain[pn] / fa;
				ib = ib * ain[pn+1] / fb;
				ic = ic * ain[pn+rowSize] / fc;
				id = id * ain[pn+rowSize+1] / fd;
				
				aout[p] += (ia+ib+ic+id);
				aout[pn] -= ia;
				aout[pn+1] -= ib;
				aout[pn+rowSize] -= ic;
				aout[pn+rowSize+1] -= id;				
			}
		}
		
		for(int p = 0; p < tsize; p++) {
			if(Math.abs(aout[p]) < 1e-8) aout[p] = 0;
		}
	}
	
	public int cell(int x, int y) {
		if(x < 0) x = 0;
		if(x >= rowSize) x = rowSize - 1;
		if(y < 0) y = 0;
		if(y >= colSize) y = colSize - 1;
		return x + y * rowSize;
	}
	
	public int ncell(int x, int y, int xn, int yn) {
		if(xn < 0) xn = 0;
		if(xn >= rowSize) xn = rowSize - 1;
		if(yn < 0) yn = 0;
		if(yn >= colSize) yn = colSize - 1;
		int pos = xn + yn * rowSize;
		if(wall[pos]) return x + y * rowSize;
		else return pos;
	}
}
