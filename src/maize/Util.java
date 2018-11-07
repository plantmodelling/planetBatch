package maize;

import javax.vecmath.Matrix3d;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.Random;

/**
 * Static class to contain in a single copy of the objects and methods useful to all other classes.
 * Avoids duplication and multiplication of objects having a potentially large memory footprint.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public final class Util {
	
	static ParameterSet parms;
		
	public static Matrix3d rTransform = new Matrix3d();
	public static final Vector3d dirHorizRef = new Vector3d(1.0,0.0,0.0);
	public static AxisAngle4d axisAngle = new AxisAngle4d(); 				// structure for rotation around an axe
	public static final Random genAleat = new Random();  					// random number generator
	public static final double epsilonDouble = 1.0e-15; 					// small double value, usefull for test
	public static final double epsilonFloat = 1.0e-7; 						// small float value, usefull for test
	public static int rootID = 0;											// root identifier
	public static Random random = new Random(System.currentTimeMillis());

	
	/**
	 * Initialize parameters
	 * @param ps
	 */
	public static void initialize(ParameterSet ps){
		parms = ps;
		log("Initialization of planetMaize.Util");
	}
    
	/**
	 * Write a string in the CrossTalk window
	 * @param string to write
	 */
	public static void log(String s) { System.out.println(s); }

	/**
	 * Log exception
	 * @param e
	 */
	public static void logException(Throwable e) {
		log(e.getClass().getName() + ": " + e.getMessage());
		StackTraceElement[] ste = e.getStackTrace();
		for (int i = 0; i < ste.length; i++) {
			StackTraceElement s = ste[i];
			log("   at " + s.getClassName() + "." + s.getMethodName() + 
					"(" + s.getFileName() + ":" + s.getLineNumber() + ")");
		}
	}
  
	/**
	 * Spin a vector around an reference from a given angle
	 * @param t: the vector to spin
	 * @param x: define the reference	
	 * @param y: define the reference
	 * @param z: define the reference
	 * @param alpha: the angle to spin
	 */
	public static void rotation3D(Vector3d t, double x, double y, double z, double alpha) {
		axisAngle.set(x, y, z, alpha);
		rTransform.set(axisAngle);
		rTransform.transform(t);
	}

	/**
	 * Spin a vector around an reference from a given angle
	 * @param t: the vector to spin
	 * @param v: the reference vector
	 * @param alpha: the angle to spin
	 */
	public static void rotation3D(Vector3d t, Vector3d v, double alpha) {
		axisAngle.set(v, alpha);
		rTransform.set(axisAngle);
		rTransform.transform(t);
	}

	
	/**
	 * Get a normalized horizontal vector which direction is random
	 * @param field: between O and 1
	 * @return the vector
	 */
	public static Vector3d randomHorizontalDirection(double field) {
		double angle = Util.genAleat.nextDouble() * 2.0 * field * Math.PI;
		return(new Vector3d(Math.cos(angle) , Math.sin(angle),0.0)); 
	}
 
  
	/**
	 * Get the volume of a sphere from its diameter
	 * @param diameter
	 * @return the volume
	 */
	public static double getSphereVolume(double diameter) { return  (Math.PI*Math.pow(diameter, 3.0) / 6.0); } 
  
	/**
	 * Get the mass of a sphere from is diameter and density
	 * @param diameter
	 * @param density
	 * @return the mass
	 */
	public static double getSphereMass(double diameter, double density) { return  (Math.PI*Math.pow(diameter,3.0)*density/6.0); } 

	/**
	 * Get the diameter of a sphere from its volume
	 * @param vol
	 * @return the diameter
	 */
	public static double getSphereDiameterFromVolume(double vol) { return  Math.pow(6.0*vol/Math.PI,1.0/3.0); } 

	
	/**
	 * Get the diameter increase of a sphere from its volume increase
	 * @param deltaVolume: the increase in volume
	 * @param diametre: the initial diameter
	 * @return the increase in diameter
	 */
	public static double getSphereDiamIncFromVolInc(double deltaVolume, double diameter) {
		return  Math.pow(6.0*(getSphereVolume(diameter)+deltaVolume)/Math.PI,1.0/3.0)-diameter;
	}

	
	/**
	 * Get the diameter increase of a sphere from its mass increase
	 * @param deltaMass: the increase in mass
	 * @param diameter: the initial diameter
	 * @param density: the density
	 * @return the diameter increase
	 */
	public static double getSphereDiamIncFromMassInc(double deltaMass, double diameter, double density) {
		return  Math.pow(6.0*(getSphereMass(diameter,density)+deltaMass)/(Math.PI*density),1.0/3.0)-diameter;
	} 
	
  	/**
   	* Get a sphere surface from its diameter
   	* @param diameter
   	* @return the surface
   	*/
  	public static double getSphereSurface(double diameter) { return  Math.PI*(diameter/2.0)*(diameter/2.0)*4;} 
  
	/**
	 * Get the surface of a cylinder (root segment, ..) from its diameter and length
	 * @param diameter of the cylinder
	 * @param length of the cylinder
	 * @return the surface
	 */
	public static double getCylinderSurface(double diameter, double length) { return  Math.PI * diameter * length; } 
  
	/**
	 * Get the volume of a cylinder (root segment, ..) from its diameter and length
	 * @param diameter of the cylinder
	 * @param length of the cylinder
	 * @return the volume
	 */
	public static double getCylinderVolume(double diameter, double length) { return  (Math.PI * diameter * diameter * length / 4.0); } 
	
  
	/**
	 * Get the length of a cylinder from its volume and diameter
	 * @param volume of the cylinder
	 * @param diametre of the cylinder
	 * @return the length
	 */
	public static double getCylinderLength(double volume, double diametre) { return  (4.0*volume/(Math.PI * diametre * diametre)); } 

	/**
	 * Get the increase in length of a cylinder from its increase in volume
	 * @param deltaVolume: the increase in volume
	 * @param diameter: the initial diameter
	 * @return the increase in length
	 */
	public static double getCylinderLengthIncFromVolInc(double deltaVolume, double diameter) { return  (4.0*deltaVolume/(Math.PI * diameter * diameter)); } 

	/**
	 * Get the increase in length of a cylinder from its increase in mass
	 * @param deltaMass: the increase in mass
	 * @param diameter: hte initial diameter
	 * @param density: the article density
	 * @return: the increase in length
	 */
	public static double getCylinderLengthIncFromMassInc(double deltaMass, double diameter, double density) {
		return  (4.0 * deltaMass / (Math.PI * diameter * diameter * density));
	}

	/**
	 * Get the increase in diameter of a cylinder from its increase in mass.
	 * Length and density are fixed
	 * @param deltaMass: the increase in mass
	 * @param length: the lenght
	 * @param diameter: the initial diameter
	 * @param density: the article density
	 * @return the increase in diameter
	 */
	public static double getCylinderDiamIncFromMassInc(double deltaMass, double length, double diameter, double density) {
		return  Math.sqrt(((Math.PI*diameter*diameter*density*length)+(4*deltaMass))/(Math.PI*length*density));
	} 

  /**
   * Get the surface of an ellipse based on its length and width
   * @param length
   * @param width
   * @return the surface
   */
	public static double getEllipseSurface(double length, double width) { return  Math.PI * length * width / 4.0f; } 

	/**
	 * Get the volume of a thick ellipse from its length, width and thickness
	 * @param length
	 * @param width
	 * @param thickness
	 * @return
	 */
	public static double getEllipseVolume(double length, double width, double thickness) { return  Math.PI * length * width * thickness/4.0f; } 

  
	/**
	 * Get the length of and ellipse from its volume, thickness and axes ratio
	 * @param volume
	 * @param axesRatio: ratio between length and width
	 * @param thickness
	 * @return the length
	 */
	public static double getEllipseLength(double volume, double axesRatio, double thickness) {
		return  Math.sqrt(4.0 * volume * axesRatio / (Math.PI * thickness));
	} 


	/**
	 * Get the increase in lenght of an ellipse from its increase in mass
	 * @param deltaMass: increase in mass
	 * @param mass: initial mass
	 * @param length: initial length
	 * @param axesRatio: ratio between length and width
	 * @param thickness
	 * @param density
	 * @return the increase in length
	 */
	public static double getEllipseLengthIncFromMassInc(double deltaMass, double mass, double length, double axesRatio, double thickness, double density) {
		return  (Math.sqrt(4.0 * axesRatio * (mass + deltaMass) / (Math.PI * thickness * density)) - length);
	} 

	/**
	 * Get the increase in thickness of an ellipse from its increase in mass
	 * @param deltaMass: increase in mass
	 * @param lenght: initial length
	 * @param axesRatio: ratio between length and width
	 * @param density
	 * @return the increase in thickness
	 */
	public static double getEllipseThickIncFromMassInc(double deltaMass, double lenght, double axesRatio, double density) {
		return  (4.0 * deltaMass * axesRatio / (Math.PI * lenght * lenght * density));
	}

  
	/**
	 * Test if a vector is vertical
	 * @param v: the vector to test
	 * @return true if vertical, false if not
	 */
	public static boolean isVertical(Vector3d v) {
		return ((Math.abs(v.x) < Util.epsilonDouble) && (Math.abs(v.y) < Util.epsilonDouble));
	} 

	
	/**
	 * Test is the given node is located inside the soil domain
	 * @param node
	 * @param soil
	 * @return true if the segment is in the soil domain
	 */
	public static boolean isOutOfSoilDomain(Point3d node, Soil soil){ 
		//if(node.x < soil.getMinX()+5 || node.x > soil.getMaxX()-5 || node.y < soil.getMinY()+5 || node.y > soil.getMaxY()-5) return true;
		return false;
	}
	
	/**
	 * Test if there it is day time (light > 0)
	 * @return true if day time
	 */
	public static boolean isDayTime(){
		return (ExogenousEnvironment.getLight(10, 10, 10, parms.maxPAR) > 0);
	}
	
	/**
	 * Get the radial angle of a vector on its "carrying vector"
	 * @param carrier: the carrying vector
	 * @param vector: the vector
	 * @return the radial angle
	 */
	public static double angleRadial(Vector3d carrier, Vector3d vector) {
		
		if (carrier.angle(vector) < Util.epsilonDouble) return 0.0; 	// radial angle is undetermined and then fixed to zero
		else {
			if (isVertical(carrier)) {
				Vector3d v1 = Util.dirHorizRef; 						// reference direction
				Vector3d v2 = new Vector3d(vector.x, vector.y, 0.0); 	// projection on the horizontale
				v2.normalize();
				Vector3d v1CrossV2 = new Vector3d();
				v1CrossV2.cross(v1,v2);
				if (v1CrossV2.dot(carrier) > 0.0) return v1.angle(v2); 
				else return (2.0 * Math.PI) - v1.angle(v2);
			}
			else {
				// Calcul of the projection of the vertical (bottom) on a plane perpendicular to the carrier
				Vector3d v1 = new Vector3d();
				v1.x = carrier.x * carrier.z / (carrier.length() * Math.sqrt((carrier.x * carrier.x) + (carrier.y * carrier.y)));
				v1.y = carrier.y * carrier.z / (carrier.length() * Math.sqrt((carrier.x * carrier.x) + (carrier.y * carrier.y)));
				v1.z = -Math.sqrt((carrier.x * carrier.x) + (carrier.y * carrier.y)) / carrier.length();
				if (carrier.dot(Util.dirHorizRef) < 0.0) v1.negate();
		
				// Calcul of the projection on a plane perpendicular to the carrier
				Vector3d v2 = new Vector3d();
				v2.x = vector.x - (carrier.x * carrier.dot(vector) / carrier.lengthSquared());
				v2.y = vector.y - (carrier.y * carrier.dot(vector) / carrier.lengthSquared());
				v2.z = vector.z - (carrier.z * carrier.dot(vector) / carrier.lengthSquared());

				Vector3d v1CrossV2 = new Vector3d();
				v1CrossV2.cross(v1,v2);
				if (v1CrossV2.dot(carrier) > 0.0) return v1.angle(v2); 
				else return (2.0 * Math.PI) - v1.angle(v2);
			}
		}
	} 

	/**
   	* Returns the value of signal in response to linear hypothesis.
   	* We assume a promoting effect of the signal (positive response, with three positive parameters).
   	* @param intercept
   	* @param slope
   	* @param signal
   	* @return the value of the response
   	*/
	public static double getLinearResponseFromPromotingSignal(double intercept, double slope, double signal) {
		if ((intercept >= 0.0) && (slope >= 0.0) && (signal >= 0.0)) return (intercept + (slope * signal)); 
		else {
			Util.log("Util.getLinearResponseFromPromotingSignal : Problem !" + " Intercept " + intercept + "/ Slope " + slope + "/ Signal " + signal);
			return 0.0001;
		}
	} 

	/**
	 * Returns the value of signal in response to linear hypothesis.
	 * We assume an inhibitory effect (or repressor) signal.
	 * (yes or no, with intercept and positive signal, negative slope).
	 * @param intercept
	 * @param slope
	 * @param signal
	 * @return the value of the response
	 */
	public static double getLinearResponseFromInhibitingSignal(double intercept, double slope, double signal) {
		if ((intercept >= 0.0) && (slope <= 0.0) && (signal >= 0.0)) return Math.max(intercept + (slope * signal),0.0); 	
		else {
			Util.log("Util.getLinearResponseFromInibitingSignal : Problem !" + " Intercept " + intercept + "/ Slope " + slope + "/ Signal " + signal);
			return 0.0;
		}
	} 

	/**
	 * Returns the value in response to the signal under the hypothesis of Michaelis dynamic (between 0 and 1)
	 * We assume a promoting effect of the signal (positive response, with both parameters positive)
	 * @param km
	 * @param signal
	 * @return the value of the response
	 */
	public static double getMichaelisResponseFromPromotingSignal(double km, double signal) {
		if ((km >= 0.0) && (signal >= 0.0)) return signal / (km + signal); 
		else {
			Util.log("Util.getMichaelisResponseFromPromotingSignal : Problem ! km : " + km + "/ signal : " + signal);
			return 0.0;
		}
	} 

	/**
	 * Returns the value in response to the signal under the hypothesis of Michaelis dynamic (between 0 and 1)
	 * We assume a promoting effect of the signal (positive response, with both parameters positive)
	 * @param km
	 * @param signal
	 * @return the value of the response
	 */
	public static double getMichaelisResponseFromInhibitingSignal(double km, double signal) {
		if ((km >= 0.0) && (signal >= 0.0)) return (1.0 - (signal / (km + signal))); 
		else {
			Util.log("Util.getMichaelisResponseFromInhibitingSignal : Problem ! km : " + km + "/ signal : " + signal);
			return 0.0;
		}
	} 
	
	/**
	 * Get the interpolated y value of a given x point, based on x and y vectors
	 * This is linear interpolation
	 * @param xVect: vector containing the x values
	 * @param yVect: vector containing the y values
	 * @param x: point from which we want the interpolated value
	 * @return the interpolated y value
	 */
	public static float getInterpolatedValue(float[] xVect, float[] yVect,float x){
		int i = 1;
		if(x < xVect[0]) return 0.0f;
		while(xVect[i] < x) i++;
		return yVect[i-1] + (x-xVect[i-1]) * ((yVect[i] - yVect[i-1]) / (xVect[i] - xVect[i-1]));
	}
  
	/**
	 * Create a new seedling, composed of a virtual article (seed), a shoot meristem and a roots meristems
	 * The seedling is positionned.
	 * @param parms
	 * @return the seedling
	 */
	public static Network createSeedling(ParameterSet parms) {

		// Create the virtual article
		Vector3d dirAx = new Vector3d(1.0, 0.0, 0.0);
		dirAx.normalize();
		Article av = new ArtVirt(new Point3d(0, 0, -2), dirAx, new Vector3d(0.0, 0.0, 1.0), new EndogenousEnvironment(0.001f,1), parms);
		Network resSemence=new Network(av);
     
		// Create an epicotyl
		Article en = new SegStem(0.6, Math.PI/2.0, 0.0, 0.0, parms.stemInitDiam, 0.01 , 1, 1, 40, parms);
		Network resEpicotyle=new Network(en);
    
		// Create a stem meristem and add it to the epicotyl
		Article mt = new MerisStem(1.0,0.0,Math.PI,0.0, parms.stemInitDiam, 1, 40, parms);
		Network resMerisTige1=new Network(mt);
		resEpicotyle.addDistalChildNetwork(resMerisTige1); 
       
		// Create an hypocotyl
		Article hc = new SegRoot(0.6, Math.PI / 2.0, Math.PI, 0.0, parms.principalRootInitDiam, 1, 10, parms, getNextRootID());
		Network resHypocotyle=new Network(hc);
    
		// Create a root meristem (primary root) and add it to the hypocotyl
		Article mr = new MerisRoot(1.0, 0.0, 0.0, Util.genAleat.nextDouble() * 2.0 * Math.PI, parms.principalRootInitDiam, 1, 10, parms, getCurrentRootID());
		Network resMerisRac=new Network(mr);
		resHypocotyle.addDistalChildNetwork(resMerisRac); 

		// Connect the epi- and hypocotyl to the virtual article
		resSemence.addDistalChildNetwork(resHypocotyle);
		resSemence.addDistalChildNetwork(resEpicotyle);
    
		// Create seminal root and connect them to the virtual article.
		// This is not completely accurate as seminal root normally emerge
		// from the epicotyl (in maize a least)
		for(int i=1; i<=parms.numberOfSeminals; i++){
    	
			// Create the segment
			double stochasticity = 1;
			if(parms.semAngStochasticity) stochasticity =  Util.genAleat.nextDouble();
			double ang1 = ((2 * Math.PI) / parms.numberOfSeminals) * (i-1) * stochasticity;
			double ang2 = ( Math.PI / i );
			Article hci = new SegRoot(0.6, ang1, ang2, 0.0, parms.seminalRootInitDiam, 1, 20, parms, getNextRootID());
			Network resHypocotylei=new Network(hci);
        	
			// Create and attach the meristem
			Article mri = new MerisRoot(1.0, 0.0 ,0.0 ,Util.genAleat.nextDouble()*2.0*Math.PI, parms.seminalRootInitDiam, 1, 20, parms, getCurrentRootID());
			Network resMerisRaci=new Network(mri);   
			resHypocotylei.addDistalChildNetwork(resMerisRaci); 

			// Connect the seminal root to the virtual article
			resSemence.addDistalChildNetwork(resHypocotylei);
		}
		// Set the position of the seedling
		resSemence.position();
		return resSemence;
	}
  
	/**
	 * Get the next root identifier
	 * @return rootID
	 */
	public static int getNextRootID(){
		parms.rootID ++;
		return parms.rootID;
	}
	
	/**
	 * Get the current root identifier
	 * @return rootID
	 */
	public static int getCurrentRootID(){ return parms.rootID; }
  
	/**
	 * Set the rootID in PlanetCTParameterSet as zero
	 */
	public static void initializeRootID(){ parms.rootID = 0; }
  
}
