package maize;

import java.io.IOException;
import java.util.ArrayList;


/**
 * Class managing the actions performed at every time step of the simulation
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Xavier Draye - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public class Simulator {
  
	// The different element of the simulation
	ParameterSet parms;
	ExogenousEnvironment envi;            
	public Network plant;            
	Time time;   
	int rewateringTime;

	ArrayList<Node> nodeList = new ArrayList<Node>();	
  	public boolean exportFirst;
	
	/**
	 * Constructor
	 */
	public Simulator() {
		this(true);
	}        
  
	/**
	 * Constructor
	 * @param interfaceNative
	 */
	public Simulator(boolean interfaceNative) {
		super();        
	}      

	/**
	 * Constructor
	 * @param interfaceNative
	 * @param parms
	 * @param soil
	 */
	public Simulator(ParameterSet p) {
		super();  
		parms = p;
		Article.initialize();
		Util.initialize(parms);
		Util.log("Welcome in Planet-Maize");
	}  
    
	/**
	 * Initialize the simulation and create the needed environment and, 
	 * if needed, the needed export tables / files
	 */
	public void initialize(){
		
		// Create the different environments
		Time.initialise(parms);   
	    envi = new ExogenousEnvironment();
	    Util.initializeRootID();
	    plant = Util.createSeedling(parms);  
	    parms.reserveC = 0;
	    rewateringTime = parms.reWateringFrequence;
	    Solver s = new Solver();
	    s.resetGrowthEff(plant);
	    
	    // Create the data tables
	    //if(parms.exportType == 1) Exports.createDataTables(parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq);
	    //if(parms.exportType == 2) Exports.createCSV(parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq);
//	    if(parms.exportType == 4) Exports.createCSVSoil(soil);
	    //Exports.exportSimulationParameters(parms, soil.parms);
	}
		
	/**
	 * Move on with the simulation
	 * @param poursuite
	 * @param nbPas
	 * @throws IOException
	 */
	public void moveOn(int nbPas) throws IOException {
		
		int i=0;

		// Initialisation of the simulation
		if ((plant==null)||(envi==null)) {
			initialize();
		}   

		while ((!Time.isOver())&&(i < nbPas)) {  
			
			//Util.log("NODE_LIST.Size : " + nodeList.size());
			Time.moveOn();
			i++;
			
			// Reset the photosynthesis value for the leaf articles
			for (Node n : nodeList){
				if(n.article instanceof SegLeaf){
					SegLeaf sl = (SegLeaf) n.article;
					sl.resetPhotosynthesis();
				}
			} 

			// Decrement the time for the next soil rewatering
			rewateringTime -= 1;
			
			// Resolutions
			for(int j = 0; j < parms.resolveIteration; j++){
  		  
				// Water potential resolution 
				// Modifications for PRTHH4
				if(parms.resolveWater){
//					Util.log("Initialization of planetMaize.Solver.solveWater()");
					Solver s = new Solver();
					if(parms.resolveWaterType == 1) s.solveWater(plant, parms);
					else s.solveWater2(plant);
				}
		
				// ABA resolution
				// Solve the ABA production 
				if(parms.ABAProduction && Time.getTime() > parms.startStress){ plant.soluteReaction(1); }
				
				// Solve the ABA fluxes
				if(parms.ABASolveMethod && Time.getTime() > parms.startStress){
					Solver s = new Solver(); 
					s.solveSoluteXylem(plant, parms, 1); 
				}
				
				// Update quantABA to the value of quantABANew
				if(parms.ABAProduction || parms.ABASolveMethod ) plant.updateSoluteValuesToNew(1);

				// Water fluxes resolution
				if(parms.resolveWater){
					// compute the water fluxes (axial and radial) based on the new water potential values
					plant.setWaterFlux();						
					// compute the water content of the soil based on the radial fluxes
//					if (Time.getTime() > parms.startStress && parmsS.updateSoil) plant.updateWaterContentExo(parms.timeStep * 3600 / parms.resolveIteration);		
				}
  	
				// Soil water resolution
//				if(j < parms.resolveIteration-1 && parmsS.updateSoil) soil.moveOn(nbPas/parms.resolveIteration);
				
				// Set the photosynthesis value for the leaf articles
				for (Node n : nodeList){
					if(n.article instanceof SegLeaf){
						SegLeaf sl = (SegLeaf)n.article;
//						sl.updatePhotosynthesis();
					}
				}  	  
				// Export data into a database or a CSV file
				if(parms.exportType > 0 && Time.getTime() > parms.startExport){
					if(parms.exportAtNoon){
						if(Time.getCurrentHour() == 12 && j == 0) exportData(parms.exportType, j);
					}
					else exportData(parms.exportType, j);
				}
				//Exports.updateExportSimulationParameters();
			}
  	  
			// Carbon resolution
			if(parms.resolveCarbon){
//				Util.log("Initialization of planetMaize.Solver.solveCarbon()");
				Solver s = new Solver(); 	
				s.solveCarbon(plant,parms); 
			}
			plant.develop();  	  
			plant.position();
			
			// Rewater the soil if needed
			if(rewateringTime <= 0){
//				soil.reWaterSoil();
				rewateringTime = parms.reWateringFrequence;
			}
			
		}   
		nodeList.clear();
		inscrit(plant, -1);
//		Util.log("--------------------------------");
	}

	/**
	 * Export the simulation data to the selected output
	 * @param type [1 = SQL, 2 = CSV, 3 = Both]
	 * @param inc: the time increment between two time steps
	 */
	public void exportData(int type, int inc){
		if(exportFirst){
			Util.log("Data export started");
			exportFirst = false;
		}
		if(type == 1){
			Exports.SQLexport(nodeList, parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq);
		}
		else if(type == 2) {
			Exports.CSVExport(nodeList, parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq, inc/parms.resolveIteration, parms);
		}
		else if(type == 3) {
			Exports.CSVExport(nodeList, parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq, (inc+1)/parms.resolveIteration, parms);
			Exports.SQLexport(nodeList, parms.exportArchitecture, parms.exportWater, parms.exportABA, parms.exportCarbon, parms.exportResEq);
		}
		else if(type == 4){
//			Exports.CSVSoilExport((inc+1)/parms.resolveIteration);
		}
	}
  
	/**
	 * Create a node list with all the article forming the plant
	 * @param r
	 * @param pere
	 */
	public void inscrit(Network r, int pere) {
		addNode(r, pere);
		int np = nodeList.size() - 1;
		if (r.childNetworkList != null) 
			for (Network rf : r.childNetworkList) 
				inscrit(rf, np);
	}
   
	/**
	 * Add node the the node list
	 * @param r
	 * @param pere
	 */
	public void addNode(Network r, int pere) {
		nodeList.add(new Node(r.baseArticle, pere));
	}

	/**
	 * Node class used for the node list containing all the articles
	 * @author Xavier Draye - Université catholique de Louvain - Earth and Life Institute (Belgium)
	 */
	class Node {
		public Article article;
		public int indicePere;
		public Node(Article a, int i) {
			article = a;
			indicePere = i;
		}
	}
}

