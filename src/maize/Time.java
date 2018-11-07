package maize;

/**
 * Class managing the time in the simulation
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 *
 */

public final class Time {

	private static double duration; 		// duration of the simulation			
	private static double timeStep; 		// time step of the simulation
	private static double time=0.0;  		// current time [hour]
	private static int day = 0;				// current day
	private static int pDay = 0;			// previous day
	private static int hour = 0;			// current hour of the day [1-24]
	private static int pHour = 0;			// previous hour of the day [1-24]
	static ParameterSet parms;		// parameter set

	
	/**
	 * Initialize the time parameters (time = 0, day = 0, ...)
	 * @param pa
	 */
	public static void initialise(ParameterSet pa) {
		time = 0;
		day = 0;
		parms = pa;
		duration = parms.simulationDuration;         
		timeStep = parms.timeStep; 
  	}

	/**
	 * Update the time parameters based on the time step value
	 */
  	public static void moveOn() {
  		time += timeStep;
  		pDay = day;
  		pHour = hour;
  		if((time-(day*24))/24 > 1) day ++;
  		hour = (int) time - (day * 24);
  	}
    
  	/**
  	 * Test weither the simulation is over
  	 * @return true is time > duration
  	 */
  	public static boolean isOver(){ return (time > duration); }

  	/**
  	 * Get the current ime of the simulation		
  	 * @return the simulation time [h]
  	 */
	public static double getTime() { return time; }

	/**
	 * Get the time step of the simulation 
	 * @return the time step [h]
	 */
	public static double getTimeStep() { return timeStep; }

	/**
	 * Get the current day of the simulation
	 * @return the current day 
	 */
	public static int getCurrentDay() { return day; }

	/**
	 * Get the current hour of the day [1-24]
	 * @return the current hour
	 */
	public static int getCurrentHour() { return hour; }

	/**
	 * Get the day of the simulation at the previous time step
	 * @return the previous day
	 */
	public static int getPreviousDay() { return pDay; }	

	/**
	 * Get the hour of the day at the previous time step
	 * @return the previous hour
	 */
	public static int getPreviousHour() { return pHour; }

}