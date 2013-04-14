package sandbox;
import processing.core.*;

public class Sand {
	PApplet parent;
	int asize;
	int w;
	int h;
	int cMat = 3;
	int[] world0;
	int[] world1;
	int[] bu;
	int[] coords;
	double[] px;
	double[] py;
	double[] vx;
	double[] vy;
	double[] fx;
	double[] fy;
	int[] mat;
	//                     noth,  wall,bouncy,  sand, water,  fire, plant,   oil,  wgas,  salt, saltw,  seed,sprout, magma, spout, stone
	boolean[] bouncy   = {false, false,  true, false, false, false, false, false, false, false, false, false, false, false, false, false};
	boolean[] solid    = {false,  true,  true,  true, false, false,  true, false, false,  true, false,  true,  true, false,  true,  true};
	boolean[] gas      = {false, false, false, false, false,  true, false, false,  true, false, false, false, false, false, false, false};
	boolean[] burnable = {false, false, false, false, false, false,  true,  true, false, false, false,  true,  true, false, false, false};
	boolean[] liquid   = {false, false, false, false,  true, false, false,  true, false, false,  true, false, false,  true, false, false};
	boolean[] still    = { true,  true, false, false, false, false,  true, false, false, false, false, false,  true, false,  true, false};
	int[] dens         = {    0,     0,     0,     0,    30,     0,     0,    10,     0,     0,    50,     0,     0,   100,     0,     0};
	byte fc = 0;
	double ax = 0;
	double ay = 10;
	//double dt = 0.1;
	double flAttr = 0.01;
	double flForce = 5000;
	double b  = (1-0.001);
	//double bn = (1-dt*0.6);
	double bouncyCol = -0.8;
	double staticCol = 0.8;
	double fluidCol = 0.1;
	
	public Sand(int wc, int hc, PApplet p) {
		parent = p;
		w = wc;
		h = hc;
		asize = w*h;
		
		world0 = new int[asize];
		world1 = new int[asize];
		coords = new int[2];
		px = new double[asize];
		py = new double[asize];
		vx = new double[asize];
		vy = new double[asize];
		fx = new double[asize];
		fy = new double[asize];
		
		int[] tmat = {0, 
				parent.color(128), 
				parent.color(255, 102, 204), 
				parent.color(238,203,128), 
				parent.color(32, 32, 255), 
				parent.color(255, 0, 0), 
				parent.color(0, 255, 0), 
				parent.color(90,60,40), 
				parent.color(230), 
				parent.color(255), 
				parent.color(180,240,255),
				parent.color(0, 200, 0), 
				parent.color(20, 200, 20), 
				parent.color(255, 100, 0), 
				parent.color(40,40,255), 
				parent.color(142)};
		
		mat = tmat;
		
		world0[200 + w*81] = 2;
		vx[200 + w*81] = 8;
		vy[200 + w*81] = 0;
		
		for(int i = 0; i < w; i++) {
			world0[i] = world1[i] = 1;
			world0[i+w*(h-1)] = world1[i+w*(h-1)] = 1;
		}
		
		for(int i = 0; i < h; i++) {
			world0[i*w] = world1[i*w] = 1;
			world0[w-1+i*w] = world1[w-1+i*w] = 1;	
		}
	}
	
	void changePhysics(double nax, double nay, double fla, double flf) {
		if(nax != 0) ax = nax;
		if(nay != 0) ay = nay;
		if(fla != 0) flAttr = fla;
		if(flf != 0) flForce = flf;
	}
	
	public int update(double dt) {
		int acc = 0;
		for(int i = 0; i < asize; i++) {
			int cMat = world0[i];
			if(cMat == 0) continue;
			if(still[cMat]) { world1[i] = cMat; continue;}
			acc++;
			double vx0;
			double vy0;
			//TODO: calc proper accelerations
			if(gas[cMat]) {
				vx0 = vx[i] = dt*(fx[i]-ax/2.0) + b*vx[i];
				vy0 = vy[i] = dt*(fy[i]-ay/2.0) + b*vy[i];
			} else {
				vx0 = vx[i] = dt*(ax+fx[i]) + vx[i];
				if(solid[world0[i+w]]) vy0 = vy[i] = dt*fy[i] + vy[i];
				else vy0 = vy[i] = dt*(ay+fy[i]) + vy[i];
			}
			if(Math.abs(vx0) < 1e-3) vx0 = vx[i] =0;
			if(Math.abs(vy0) < 1e-3) vy0 = vy[i] =0;
			double px0 = px[i] += dt*vx[i];
			double py0 = py[i] += dt*vy[i];
			int dx = (int) px[i];
			int dy = (int) py[i];
			int pn = i + dx + dy*w;
			//PApplet.println(vx0 + " " +vy0 + " " +px0 + " " +py0 + " " +dx + " " +dy+ " " +fx[i] + " " +fy[i]);
			if(pn == i) {
				world1[i] = cMat;
				continue;
			}
			
			if(Math.abs(dx) > 1 || Math.abs(dy) > 1) {
				//moving more than one pixel
				//check if we hit something on the way
				int x0 = i % w;
				int x1 = x0 + dx;
				int y0 = i / w;
				int y1 = y0 + dy;
				pathCollisions(x0, y0, x1, y1);
				//println(x0 + ", " + y0 + ", " + x1 + ", " + y1 + ", " + coords[0] + ", " + coords[1]);
				int npn = coords[0] + w*coords[1];
				if(pn != npn) {
					//we hit something
					if(bouncy[cMat]) {
						//TODO: boundaries on pns!
						if(world0[npn+1] != 0 || world1[npn+1] != 0 || world0[npn-1] != 0 || world1[npn-1] != 0)
							vx0 *= bouncyCol;
						if(world0[npn+w] != 0 || world1[npn+w] != 0 || world0[npn-w] != 0 || world1[npn-w] != 0)
							vy0 *= bouncyCol;
						
					}
					else {
						vx0 *= staticCol;
						vy0 *= staticCol;
					}
				}
				pn = npn;
				if(pn > asize) pn = pn - asize;
				if(pn < 0) pn = pn + asize;
			} else {
				if(pn > asize) pn = pn - asize;
				if(pn < 0) pn = pn + asize;
				//just check next pixel
				/********************************************
				 * There is still no mass conservation
				 * because if a pixel stays still where someone
				 * already is in world1 then one of them disappears
				 *******************************************/
				if(world0[pn] != 0 || world1[pn] != 0) {
					if((solid[cMat] || liquid[cMat]) && gas[world0[pn]]){
						world0[i]= world1[i] = world0[pn];
						world0[pn] = 0;//empty
						vx0 *= staticCol;
						vy0 *= staticCol;
					} else if((solid[cMat] || liquid[cMat]) && gas[world1[pn]]) {
						world0[i]= world1[i] = world1[pn];
						world0[pn] = 0;//empty
						vx0 *= staticCol;
						vy0 *= staticCol;
					} else if(solid[cMat] && liquid[world0[pn]]){
						//sink
						world0[i]= world1[i] = world0[pn];
						world0[pn] = 0;//empty
						vx0 *= fluidCol;
						vy0 *= fluidCol;
					} else if(solid[cMat] && liquid[world1[pn]]){
						//sink
						world0[i]= world1[i] = world1[pn];
						world0[pn] = 0;//empty
						vx0 *= fluidCol;
						vy0 *= fluidCol;
					} else if(liquid[cMat] && liquid[world0[pn]] && (dens[cMat] > dens[world0[pn]])) {
						//sink
						world0[i]= world1[i] = world0[pn];
						world0[pn] = 0;//empty
						vx0 *= fluidCol;
						vy0 *= fluidCol;
					} else if(liquid[cMat] && liquid[world1[pn]] && (dens[cMat] > dens[world1[pn]])) {
						//sink
						world0[i]= world1[i] = world1[pn];
						world0[pn] = 0;//empty
						vx0 *= fluidCol;
						vy0 *= fluidCol;
					} else {
						pn = i;
						if(bouncy[cMat]) {
							//TODO: boundaries on pns!
							if(world0[pn+1] != 0 || world1[pn+1] != 0 || world0[pn-1] != 0 || world1[pn-1] != 0)
								vx0 *= bouncyCol;
							if(world0[pn+w] != 0 || world1[pn+w] != 0 || world0[pn-w] != 0 || world1[pn-w] != 0)
								vy0 *= bouncyCol;
						}
						else {
							vx0 *= staticCol;
							vy0 *= staticCol;
						}
					}
				}
			}
			
			world1[pn] = cMat;
			px[pn] = px0 - dx;
			py[pn] = py0 - dy;
			vx[pn] = vx0;
			vy[pn] = vy0;
		}
		
		for(int i = 0; i < asize; i++) {
			int cMat = world1[i];
			if(cMat == 0 || cMat == 1) continue;
			if(liquid[cMat]) {
				//vx[i] += dt*2*(Math.random()*2 - 1);
				if(world1[i+1] != 0) vx[i] -= dt*8;
				if(world1[i-1] != 0) vx[i] += dt*8;
			}
			if(gas[cMat]) {
				vx[i] += dt*5*(Math.random()*2 - 1);
			}
			if(cMat == 5) { //fire
				if(Math.random() < 0.04) world1[i] = 0;
				if(burnable[world1[i+1]] && Math.random() < 0.1) world1[i+1] = 5;
				if(burnable[world1[i-1]] && Math.random() < 0.1) world1[i-1] = 5;
				if(burnable[world1[i+w]] && Math.random() < 0.1) world1[i+w] = 5;
				if(burnable[world1[i-w]] && Math.random() < 0.1) world1[i-w] = 5;
				if(world1[i+1] == 4 && Math.random() < 0.8) { world1[i+1] = 8; world1[i] = 0; }
				if(world1[i-1] == 4 && Math.random() < 0.8) { world1[i-1] = 8; world1[i] = 0; }
				if(world1[i+w] == 4 && Math.random() < 0.8) { world1[i+w] = 8; world1[i] = 0; }
				if(world1[i-w] == 4 && Math.random() < 0.8) { world1[i-w] = 8; world1[i] = 0; }
			}
			if(cMat == 6) { //plant
				if(world1[i+1] == 4 && Math.random() < 0.02) { world1[i+1] = 0; world1[i] = 12; }
				if(world1[i-1] == 4 && Math.random() < 0.02) { world1[i-1] = 0; world1[i] = 12; }
				if(world1[i+w] == 4 && Math.random() < 0.02) { world1[i+w] = 0; world1[i] = 12; }
				if(world1[i-w] == 4 && Math.random() < 0.02) { world1[i-w] = 0; world1[i] = 12; }
			}
			if(cMat == 8) { //wgas
				if(solid[world1[i-w]] && Math.random() < 0.001) world1[i] = 4;
			}
			if(cMat == 9) { //salt
				if(world1[i+1] == 4) { world1[i+1] = 10; world1[i] = 0; }
				if(world1[i-1] == 4) { world1[i-1] = 10; world1[i] = 0; }
				if(world1[i+w] == 4) { world1[i+w] = 10; world1[i] = 0; }
				if(world1[i-w] == 4) { world1[i-w] = 10; world1[i] = 0; }
			}
			if(cMat == 10) { //swater
				if(world1[i+1] == 6 && Math.random() < 0.05) world1[i+1] = 0;
				if(world1[i-1] == 6 && Math.random() < 0.05) world1[i-1] = 0;
				if(world1[i+w] == 6 && Math.random() < 0.05) world1[i+w] = 0;
				if(world1[i-w] == 6 && Math.random() < 0.05) world1[i-w] = 0;
			}
			if(cMat == 11) { //seed
				if(world1[i+1] == 4 && Math.random() < 0.03) { world1[i+1] = 11; world1[i] = 6; }
				if(world1[i-1] == 4 && Math.random() < 0.03) { world1[i-1] = 11; world1[i] = 6; }
				if(world1[i+w] == 4 && Math.random() < 0.03) { world1[i+w] = 11; world1[i] = 6; }
				if(world1[i-w] == 4 && Math.random() < 0.03) { world1[i-w] = 11; world1[i] = 6; }
				if(world1[i+1] == 6 && world1[i-1] == 6 && world1[i+w] == 6 &&
						world1[i-w] == 6) world1[i] = 6;
			}
			if(cMat == 12) { //sprout
				double r1 = Math.random();
				if(world1[i+1] == 6 && r1 < 0.3) { world1[i+1] = 12; world1[i] = 6; }
				if(world1[i-1] == 6 && r1 > 0.3 && r1 < 0.6) { world1[i-1] = 12; world1[i] = 6; }
				if(world1[i-w] == 6 && r1 > 0.6 && r1 < 0.9) { world1[i-w] = 12; world1[i] = 6; }
				if(!solid[world1[i-w]] && Math.random() < 0.02){ world1[i-w] = 6; world1[i] = 6; }
				if(!solid[world1[i-1]] && Math.random() < 0.01){ world1[i-1] = 6; world1[i] = 6; }
				if(!solid[world1[i+1]] && Math.random() < 0.01){ world1[i+1] = 6; world1[i] = 6; }
				/*if(world1[i+1] != 6 && world1[i-1] != 6 && world1[i+w] != 6 &&
						world1[i-w] != 6 && Math.random() < 0.08) world1[i] = 6;*/
			}
			if(cMat == 13) { //magma
				if(burnable[world1[i+1]] && Math.random() < 0.1) world1[i+1] = 5;
				if(burnable[world1[i-1]] && Math.random() < 0.1) world1[i-1] = 5;
				if(burnable[world1[i+w]] && Math.random() < 0.1) world1[i+w] = 5;
				if(burnable[world1[i-w]] && Math.random() < 0.1) world1[i-w] = 5;
				if(world1[i+1] == 4 && Math.random() < 0.8) { world1[i+1] = 8; world1[i] = 15; }
				if(world1[i-1] == 4 && Math.random() < 0.8) { world1[i-1] = 8; world1[i] = 15; }
				if(world1[i+w] == 4 && Math.random() < 0.8) { world1[i+w] = 8; world1[i] = 15; }
				if(world1[i-w] == 4 && Math.random() < 0.8) { world1[i-w] = 8; world1[i] = 15; }
				if(world1[i+1] == 3 && Math.random() < 0.02) { world1[i+1] = 13; }
				if(world1[i-1] == 3 && Math.random() < 0.02) { world1[i-1] = 13; }
				if(world1[i+w] == 3 && Math.random() < 0.02) { world1[i+w] = 13; }
				if(world1[i-w] == 3 && Math.random() < 0.02) { world1[i-w] = 13; }
			}
			if(cMat == 14) { //spout
				if(world1[i+1] == 0 && Math.random() < 0.05) { world1[i+1] = 4; }
				if(world1[i-1] == 0 && Math.random() < 0.05) { world1[i-1] = 4; }
				if(world1[i+w] == 0 && Math.random() < 0.05) { world1[i+w] = 4; }
				if(world1[i-w] == 0 && Math.random() < 0.05) { world1[i-w] = 4; }
			}
		}
		
		bu = world0;
		world0 = world1;
		world1 = bu;
		
		return acc;
	}
	
	public void fly(int cw, int ch, FluidCom fluid) {
		float[] vxf = fluid.getVx0();
		float[] vyf = fluid.getVy0();
		int[] numParticles = fluid.getNumParticles();
		int cellw = Math.round(w / (cw - 2));
		int cellh = Math.round(h / (ch - 2));
		for (int i = 0; i < asize; i++) {
			int cMat = world0[i];
			fx[i] = fy[i] = 0;
			if(cMat == 0) continue;
			int pvx = (int) (i % w / cellw) + 1;
			int pvy = (int) (i / h / cellh) + 1;
			float vvx = fluid.getMACCellX(pvx,pvy,vxf);
			float vvy = fluid.getMACCellY(pvx,pvy,vyf);
			int p = pvx + pvy * cw;
			if(still[cMat]) numParticles[p] += 1;
			else {
				//PApplet.println(i + " " + pvx + " " + pvy);
				double vxn = flForce*vvx - vx[i];
				double vyn = flForce*vvy - vy[i];
				fx[i] = flAttr*vxn;
				fy[i] = flAttr*vyn;
			}
			if(cMat == 5) {
				fluid.addTemp(i % w, i / h, 0.01f);
			}
		}
	}
	
	public void draw() {
		parent.loadPixels();
		for(int i = 0; i < asize; i++) {
			int t = mat[world0[i]];
			if(t != 0) parent.pixels[i] = t;
			world1[i] = 0;
		}
		parent.updatePixels();
	}
	
	public void addMat(int mouseX, int mouseY, int mouseXp, int mouseYp) {
		drawLine(mouseXp, mouseYp, mouseX, mouseY, cMat, world0);
		drawLineDouble(mouseXp, mouseYp, mouseX, mouseY, cMat, px, 0);
		drawLineDouble(mouseXp, mouseYp, mouseX, mouseY, cMat, py, 0);
		drawLineDouble(mouseXp, mouseYp, mouseX, mouseY, cMat, vx, 1);
		drawLineDouble(mouseXp, mouseYp, mouseX, mouseY, cMat, vy, 0);
		mouseXp = mouseX;
		mouseYp = mouseY;
	}
	
	public void setcMat(int key) {
		switch(key) {
			case 49:
				cMat = 1;//wall
				break;
			case 98:
				cMat = 2;//bouncy
				break;
			case 115:
				cMat = 3;//sand
				break;
			case 119:
				cMat = 4;//water
				break;
			case 102:
				cMat = 5;//fire
				break;
			case 112:
				cMat = 11;//seed
				break;
			case 111:
				cMat = 7;//oil
				break;
			case 103:
				cMat = 8;//gas
				break;
			case 100:
				cMat = 9;//salt
				break;
			case 109:
				cMat = 13;//magma
				break;
			case 101:
				cMat = 14;//spout
				break;
		}
	}
	
	private void pathCollisions(int x0, int y0, int x1, int y1) {
		int dx = Math.abs(x1-x0);
		int dy = Math.abs(y1-y0);
		int sx = -1;
		int sy = -1;
		if (x0 < x1) sx = 1;
		if (y0 < y1) sy = 1;
		int err = dx-dy;
		int x0p = x0;
		int y0p = y0;
		 
		while(true) {
			if (x0 == x1 && y0 == y1) {
				//finished the line without hitting something
				coords[0] = x1;
				coords[1] = y1;
				break;
			}
			int e2 = 2*err;
			if (e2 > -dy) { 
				err = err - dy;
				x0p = x0;
				x0 = x0 + sx;
			}
			if (e2 < dx) { 
				err = err + dx;
				y0p = y0;
				y0 = y0 + sy;
			}
			if(world0[x0 + y0*w] != 0 || world1[x0 + y0*w] != 0) {
				//something is in the way
				//so return previous open point
				coords[0] = x0p;
				coords[1] = y0p;
				break;
			}
		}
	}
	
	private void drawLine(int x0, int y0, int x1, int y1, int c, int[] t) {
		int dx = Math.abs(x1-x0);
		int dy = Math.abs(y1-y0);
		int sx = -1;
		int sy = -1;
		if (x0 < x1) sx = 1;
		if (y0 < y1) sy = 1;
		int err = dx-dy;
		 
		while(true) {
			int np = x0 + w*y0;
			if(np > 0 && np < asize && t[np] == 0) {
				t[np] = c; //setPixel(x0,y0);
			}
			if (x0 == x1 && y0 == y1) break;
			int e2 = 2*err;
			if (e2 > -dy) { 
				err = err - dy;
				x0 = x0 + sx;
			}
			if (e2 < dx) { 
				err = err + dx;
				y0 = y0 + sy;
			}
		}
	}
	
	private void drawLineDouble(int x0, int y0, int x1, int y1, int c, double[] t, double r) {
		int dx = Math.abs(x1-x0);
		int dy = Math.abs(y1-y0);
		int sx = -1;
		int sy = -1;
		if (x0 < x1) sx = 1;
		if (y0 < y1) sy = 1;
		int err = dx-dy;
		 
		while(true) {
			int np = x0 + w*y0;
			if(np > 0 && np < asize) {
				t[np] = Math.random() * 2 * r - r; //setPixel(x0,y0);
			}
			if (x0 == x1 && y0 == y1) break;
			int e2 = 2*err;
			if (e2 > -dy) { 
				err = err - dy;
				x0 = x0 + sx;
			}
			if (e2 < dx) { 
				err = err + dx;
				y0 = y0 + sy;
			}
		}
	}
	
	public void reset() {
		for(int i = 0; i < asize; i++) {
			world0[i] = world1[i] = 0;
			px[i] = py[i] = vy[i] = vx[i] = 0;
		}
		for(int i = 0; i < w; i++) {
			world0[i] = world1[i] = 1;
			world0[i+w*(h-1)] = world1[i+w*(h-1)] = 1;
		}
		
		for(int i = 0; i < h; i++) {
			world0[i*w] = world1[i*w] = 1;
			world0[w-1+i*w] = world1[w-1+i*w] = 1;	
		}
	}
}
