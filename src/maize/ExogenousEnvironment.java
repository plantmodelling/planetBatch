package maize;

import javax.vecmath.Point3d;

/**
 * Class pointing to the right environment, either soil or atmosphere and returning
 * environmental values.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 *
 */

public class ExogenousEnvironment {
	static Air air;				// the atmosphere
	static Soil soil;		// the soil

	
	/**
	 * Contructor
	 * @param SoilCTParameterSet
	 */
	public ExogenousEnvironment() {
		super();
		air = new Air();
		soil = new Soil();
	}
  
	/**
	 * Return the sum of temperature since the beginnning of the simulation
	 * @return sum of temperature [°C]
	 */
	public double getDeltaST(float maxT){
		double cT; 
		double pT;
		double sT1 = 0;
		double sT2 = 0;
	  
		for(int i = 1 ; i <= Time.getCurrentHour() ; i ++){
			cT = getTemperature(i ,null, null, maxT); 		// current temperature
			pT = getTemperature(i-1 ,null, null, maxT); 	// previous temperature
			sT1 += (cT+pT)/2;  
		}
		for(int i = 1 ; i <= 24 ; i ++){
			cT = getTemperature(i ,null, null, maxT); 		// current temperature
			pT = getTemperature(i-1 ,null, null, maxT); 	// previous temperature
			sT2 += (cT+pT)/2;  
		}	  
		return sT1 + (sT2 * (Time.getCurrentDay()));
	}
  
	/**
	 * Is the given position in the soil? 
	 * @param position of the article
	 * @return true/false
	 */
	public static boolean isInSoil(Point3d position) {return position.z < 0.0; }
  
	/**
	 * Is the given position in the air? 
	 * @param position of the article
	 * @return true/false
	 */
	public static boolean isInAir(Point3d position){return position.z >= 0.0;}
  
	/**
	 * Is the given position in the soil, between a given position and the surface? 
	 * @param position of the article
	 * @param depth
	 * @return true/false
	 */
	public static boolean isInSoilNearSurface(Point3d position, double depth) {
		if ((position.z < 0.0) && (position.z > -depth)) return true;
		else return false; 
	}
  
	/**
	 * Is the given position in the air, between a given position and the surface? 
	 * @param position of the article
	 * @param height
	 * @return true/false
	 */
	public static boolean isInAirNearSurface(Point3d position, double height) {
		if ((position.z > 0.0) && (position.z < height)) return true;
		else return false; 
	}
  
	/**
	 * Get the mechanical constrain of the environment based on a given position.
	 * If in the air, will return 0. 
	 * @param position of the article
	 * @param sti of the article
	 * @return mechanical constrain [0-1]
	 */
	public static double getMechanicalConstrain(Point3d position, SpaceTimeInfo sti) {
		if (isInSoil(position)) return 0;// TODO soil.getDouble(sti);
		else return 0; 
	}
  
	/**
	 * Get the temperature of the environment based on a given position.
	 * @param position of the article
	 * @param sti of the article
	 * @return temperature [°C]
	 */
	public static double getTemperature(Point3d position, SpaceTimeInfo sti) {
		if(isInSoil(position)) return Soil.getTemperature(position);
		else return Air.getTemperature(position); 
	}
  
	/**
	 * Get the temperature of the environment based on a given position and time.
	 * @param time
	 * @param position of the article
	 * @param sti of the article
	 * @return temperature [°C]
	 */	
	public static double getTemperature(int time, Point3d position, SpaceTimeInfo sti, float maxT) {
		if(isInSoil(position)) return Soil.getTemperature(time, position);
		else return Air.getTemperature(time, position, maxT); 
	}

	/**
	 * Get the VPD of the environment based on a given position.
	 * If in the soil, will return 0. 
	 * @param position of the article
	 * @return VPD [MPa]
	 */	
	public static double getVPD(Point3d position) {
		if(isInSoil(position)) return 0;
		else return Air.getVPD(position); 
	}

	/**
	 * Get the water potential of the environment based on a given position.
	 * For the needs of some simulations, some of the following method gives arbitrary values
	 * of water potential, not depending on the soil model. the type of stress is 
	 * defined by the user in CrossTalk
	 * @param position of the article
	 * @param stressType: the type of stress (homogenous, split, ...)
	 * @param stressStart: the beginning of the stress
	 * @return water potential [MPa]
	 */	
	public static double getWaterPotExo(Point3d position, int stressType, int stressStart, int rewatering) {
		  
		if(isInSoil(position)){
	  		  
			int maxT = 50;
			float slope = 0.015f;
			double pseudoWPInit = -0.15;
			double pseudoWP = pseudoWPInit;
			double minPot = -1 * slope * maxT + pseudoWPInit;
			  			  
			if (Time.getTime() > stressStart){
				double t = Time.getTime() - stressStart;
				
				switch (stressType){				
				case 0: 
					break;
			  			  
				case 1: // Soil drying homogeniously as a function of time then rewatering
					double minPot2 =  -1 * slope * 0.5 * maxT + pseudoWPInit;
					if (t < maxT) pseudoWP = -1 * slope * 0.5 * t + pseudoWPInit;
					else if (t >= maxT && t < (2*maxT)) pseudoWP = slope * 0.5 * (t-maxT) + minPot2;
					else pseudoWP = pseudoWPInit;
					break;
			  			  
				case 2: // Soil drying homogeniously as a function of time then rewatering - quick
					if (t < maxT/10) pseudoWP = -1 * slope * 10 * t + pseudoWPInit;
					else if (t >= maxT/10 && t < 100 - (maxT / 10)) pseudoWP = minPot;
					else if (t >= 100 - (maxT / 10) && t < 100) pseudoWP = slope * 10 * (t-(100 - (maxT / 10))) + minPot;
					else pseudoWP = pseudoWPInit;
					break;
			  			  
				case 3: // Vertical Partial rootzone drying and rewatering
					if (position.x > 0)  pseudoWP = pseudoWPInit;   												// One part of the soil provides enough water
					else if (position.x <= 0 && t < maxT) pseudoWP = -1 * slope * t + pseudoWPInit; 				// Second part of the soil is drying
					else if (position.x <= 0 && t >= maxT && t < (2*maxT)) pseudoWP = slope * (t-maxT) + minPot;
					else pseudoWP = pseudoWPInit;
					break;
						  
				case 4: // Vertical Partial rootzone drying quick drying then rewatering
					if (position.x > 0)  pseudoWP = pseudoWPInit;   												// One part of the soil provides enough water
					else if (position.x <= 0 && t < maxT / 10) pseudoWP = -1 * slope * 10 * t + pseudoWPInit;		// Second part of the soil is drying
					else if (position.x <= 0 && t >= maxT / 10 && t < 100-(maxT/10)) pseudoWP = minPot;
					else if (position.x <= 0 && t > 100-(maxT/10) && t < 100) pseudoWP = slope * 10 * (t-(100 - (maxT / 10))) + minPot;
					else pseudoWP = pseudoWPInit;
					break;
						  
				case 5: // Horizontal Partial rootzone drying top and rewatering
					if (position.z < -10)  pseudoWP = pseudoWPInit;   												// One part of the soil provides enough water
					else if (position.z >= -10 && t < maxT) pseudoWP = -1 * slope * t + pseudoWPInit; 				// Second part of the soil is drying
					else if (position.z >= -10 && t >= maxT && t < (2*maxT)) pseudoWP = slope * (t-maxT) + minPot;
					else pseudoWP = pseudoWPInit;
					break;
						  
				case 6: // Horizontal Partial rootzone drying bottom and rewatering
					if (position.z > -10)  pseudoWP = pseudoWPInit;   											// One part of the soil provides enough water
					else if (position.z <= -10 && t < maxT) pseudoWP = -1 * slope * t + pseudoWPInit; 			// Second part of the soil is drying
					else if (position.z <= -10 && t >= maxT && t < (2*maxT)) pseudoWP = slope * (t-maxT) + minPot;
					else pseudoWP = pseudoWPInit;
					break;
						  
				case 7: // Vertical Partial rootzone drying quick drying
					if (position.x > 0)  pseudoWP = pseudoWPInit;   					// One part of the soil provides enough water
					else pseudoWP = -1 * slope * t + pseudoWPInit; 				// Second part of the soil is drying
					break;
						  
				case 8: // Horizontal Partial rootzone drying quick drying TOP
					if (position.z > -10)  pseudoWP = pseudoWPInit;   				// One part of the soil provides egnough water
					else pseudoWP = -1 * slope * t + pseudoWPInit; 					// Second part of the soil is drying
					break;
						  
				case 9: // Horizontal Partial rootzone drying quick drying BOTTOM
					if (position.z <= -10)  pseudoWP = pseudoWPInit;   				// One part of the soil provides egnough water
					else pseudoWP = -1 * slope * t + pseudoWPInit; 					// Second part of the soil is drying
					break;
						  
				case 10: // Homogenous drying
					pseudoWP = -1 * slope * 0.5 * t + pseudoWPInit; 
					break;	  					  	
				
				case 11: // vertical gradient
					double b = - pseudoWPInit;
					double d = - position.z;
					double a = b * (t/15);
					pseudoWP = -(a + -((a - b) / 30) * d);
					if(pseudoWP > pseudoWPInit) pseudoWP = pseudoWPInit;
					if(d > 30) pseudoWP = pseudoWPInit;
					if(Time.getTime() >= rewatering) pseudoWP = pseudoWPInit;
//					System.out.println(Time.getTime());
//					System.out.println(pseudoWP);
//					System.out.println("-------");
					break;
				}
				
			}
			return pseudoWP;
		}
		else return Air.getWaterPot(position); 
	}

	/**
	 * Get the light of the environment based on a given position.
	 * If in the soil, will return 0. 
	 * @param position of the article
	 * @param sti of the article
	 * @return light [µmol/m2/sec]
	 */
	public static double getLight(Point3d position, float maxPAR) {
		if(isInSoil(position)) return 0;
		return Air.getLight(position, maxPAR); 
	}
	
	/**
	 * Get the light of the environment based on a given position.
	 * If in the soil, will return 0. 
	 * @param x
	 * @param y
	 * @param z
	 * @return light [µmol/m2/sec]
	 */
	public static double getLight(double x, double y, double z, float maxPAR) {
		return getLight(new Point3d(x, y, z), maxPAR); 
	}
}
