/* Created on May 30, 2005 */
package maize;

/** @author Xavier Draye - Universit� catholique de Louvain (Belgium) */
public class Planet {

   /*
    * run the simulation in batch mode
    * args[0] = length of the simulation
    */
   public static void main(String[] args) {
	     
	// Initialise the parameter set
	 ParameterSet ps = new ParameterSet(args[0]);        	 
	 ps.updateParameterValues();
	 
	 Simulator plant = new Simulator(ps);
	 plant.initialize();
	  
	 //for (int time = 0; time <= Integer.valueOf(args[0]); time ++) {
	 for (int time = 0; time <= ps.simulationDuration; time ++) {
		 try {plant.moveOn(1);}
		 catch(Exception e) { Util.logException(e); }
	 }    
	 System.out.println("--------------------------------------");
	 System.out.println("Simulation done. \n Final number of nodes: "+plant.nodeList.size()+" \n Simulation time: "+ps.simulationDuration);
   }
}

