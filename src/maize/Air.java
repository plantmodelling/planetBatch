package maize;

import javax.vecmath.Point3d;

/**
 * Class describing the atmosphere. Has to be complete to follow the same structure as the soil (SpaceTimeInfo) TODO
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain  - Earth and Life Institute (Belgium)
 *
 */

public class Air{

	/**
	 * Constructor
	 */
	public Air() {
	}
	
	/**
	 * Return the temperature as a function of space and time
	 * From personnal data (Louvain-la-Neuve, Belgium, greenhouses: 50°39'58.21"N - 4°37'10.44"E elev 142m)
	 * @param position
	 * @return temperature [°C]
	 */
	public static double getTemperature(Point3d position) {
	
		double temp;
		int time = Time.getCurrentHour();  
		
		if(time < 5) temp = 20;
		else temp = (-0.3*time*time)+(8*time)-25;
		if(temp < 20) temp = 20;
		  
		return temp;
	}
	
	/**
	 * Return the temperature as a function of space and time
	 * From personnal data (Louvain-la-Neuve, Belgium, greenhouses: 50°39'58.21"N - 4°37'10.44"E elev 142m)
	 * @param time
	 * @param position
	 * @return temperature [°C]
	 */
	public static double getTemperature(int time, Point3d position, float maxT) {
	
		double temp;
		
		if(time < 5) temp = 20;
		else temp = (-0.3*time*time)+(8*time)-(55-maxT);
		if(temp < 20) temp = 20;
		  
		return temp;
	}
	
	/**
	 * Return the humidity as a function of space and time
	 * From personnal data (Louvain-la-Neuve, Belgium, greenhouses: 50°39'58.21"N - 4°37'10.44"E elev 142m)
	 * @param position
	 * @return humidity [%]
	 */	  
	public static double getHumidity(Point3d position) {
	
		double hum;
		int time = Time.getCurrentHour();  
		
		if(time < 5) hum = 50;
		else hum = (0.35*time*time)-(10*time)+94;
		
		if(hum > 50) hum = 50;
		
		return hum;
	}

	/**
	 * Return the VPD as a function of space and time
	 * From personnal data (Louvain-la-Neuve, Belgium, greenhouses: 50°39'58.21"N - 4°37'10.44"E elev 142m)
	 * @param position
	 * @return VPD [MPa]
	 */	
	public static double getVPD(Point3d position) {
	
		double vpd;
		int time = Time.getCurrentHour();  
	
		if(time > 6 && time < 20) vpd = (-0.1 * time * time)+(2.7 * time) - 13;
		else vpd = 0;
			
		return vpd;
	}

	/**
	 * Return the light as a function of space and time
	 * From personnal data (Louvain-la-Neuve, Belgium, greenhouses: 50°39'58.21"N - 4°37'10.44"E elev 142m)
	 * @param position
	 * @return light [µmol/m2/sec]
	 */		
	public static double getLight(Point3d position, float maxPAR) {
	
		double light;
		int time = Time.getCurrentHour();  
			
		if (time > 5 && time < 20) light = (-7 * time * time)+(184 * time) - (1100 - maxPAR);
		else light = 0;
		
	    return light; 
	}
	 
	/**
	 * Return the water potential as a function of space and time
	 * From Girardin 98 (fig.148 p.134)
	 * @param position
	 * @return water potential [MPa]
	 */		 
	public static double getWaterPot(Point3d position) {
		  
		   double wp = 0;
			int time = Time.getCurrentHour();  
		   
		   if (time > 6 && time < 20) wp = (1.03 * time * time) - (27.6 * time) + 93;
		   else wp = 0;
		  
		   return wp;  
	}

}
