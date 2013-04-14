package sandbox;
import processing.core.*;

public class MainBox extends PApplet {
	
	private static final long serialVersionUID = -6185508289654966583L;
	int width = 400;
	int height = 400;
	int hcells = 80;
	int vcells = 80;
	int pmx;
	int pmy;
	int fmx;
	int fmy;
	int pnum = 0;
	boolean go = true;
	float angle = (float) 5.2;
	FluidCom fluid;
	Sand sand;
	
	public void setup() {
		size(width,height,P2D);
		
		fluid = new FluidCom(hcells, vcells, width, height, this, 500.0f);
		sand = new Sand(width, height, this);
	}
	
	public void draw() {
		
		if(go) {	
			float a = 0.6f/ (float) frameRate;
			fluid.update(a);
			fluid.updateWalls();
			sand.fly(hcells, vcells, fluid);
			pnum = sand.update(a);
		}
		
		fluid.draw();
		sand.draw();
		
		stroke(255,0,0);
		line((float)mouseX, (float)mouseY, (float)mouseX + 30*cos(angle), (float)mouseY + 30*sin(angle));
		stroke(0);
		
		fill(255);
		text((int)frameRate + " " + pnum,20,20);
		
		if (mousePressed) {
			//if(mouseButton == RIGHT) fluid.addDens(mouseX, mouseY);
			if(mouseButton == RIGHT) fluid.addVel(mouseX, mouseY, angle);
			if(mouseButton ==  LEFT) sand.addMat(mouseX, mouseY, pmx, pmy);
		}
		
		pmx = mouseX;
		pmy = mouseY;
	}
	
	public static void main(String args[]) {
		PApplet.main(new String[] { "--present", "MainBox" });
	}
	
	public void mousePressed() {
		fmx = mouseX;
		fmy = mouseY;
	}
	
	public void mouseMoved() {
		fmx = mouseX;
		fmy = mouseY;
	}
	
	public void keyPressed() {
		println("key event detected with num: " + (int)key);
		switch(key) {
		case 32:
			if(go) go = false;
			else go = true;
			break;
		case 110:
			angle += 0.2;
			println(angle);
			break;
		case 109:
			angle -= 0.2;
			println(angle);
			break;
		case 114:
			sand.reset();
			fluid.reset();
			break;
		case 118://v
			fluid.toggleVelocities();
			break;
		default:
			sand.setcMat(key);
			break;
		}
	}
}
