package maize;

/**
 * This class represent a stem meristem, responsible for the creation of leaves and adventious roots
 * (depending on the stem node)
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium) *
 */

public class MerisStem extends Article {

	// Architecture
	public int nodeLevel = 2;				// The node level on the plant. Start at 2 because the first node is build by the seed.
	public int leafNum = 1;					// The leaf number
	public boolean firstLeaf = true;		// Flag for the node bearing the first leaf
	public static double dr0 = 0.014;		// The minimal diameter for root initiation
	public double angPhylloAdvRoot;			// Angle for the production of adventious roots
	public double nbWaitingPlasto=0.0; 		// Number of plastochon waiting
	public double angPhyllo = 0.0; 			// Phyllotaxic angle  
	public int maxAdvRoot = 0;				// Maximal number of adventive roots
	
	// Growth and developemnt, modified from Drouet and Pagès 2003
	public double plasto;
	public double arp = 0.0349;
	public double brp = 0.001;
	public double crp = 4.208e-6;				
	public double drp = -5.208e-9; 			
  
	// ABA
	public float ABAHalfLife = 0.7f;		// Half life of ABA [h] (Ren et al. 2007)
	
	// Carbon
	private double dryMassRatio = 0.01f;							// The dry mass ration of the article [gC/cm3]
	public double cReserve = 0;										// The reserve in C for the production of new organs
	public double stemEff = 1, leafEff = 1, rootEff = 1;			// The growth efficiency coefficients for stem, leaf and roots
	public double stemDemand = 0, leafDemand = 0, rootDemand = 0;	// The growth demand for stem, leaf and roots
	public double growthDemand = 0;									// The total growth demand 
	
	/**
	 * Constructor
	 */
	public MerisStem() {
		super();
	}

	/**
	 * Constructor
	 * @param posRel
	 * @param angIns
	 * @param angRad
	 * @param angSpin
	 * @param diam
	 * @param ordre
	 * @param type
	 * @param parms
	 */
	public MerisStem(double posRel, double angIns, double angRad, double angSpin,
			double diam, int ordre, int type, ParameterSet parms) {
		super(posRel, angIns, angRad, angSpin, ordre, type, parms);
		diameter=diam;
		angPhylloAdvRoot = 0;
		setParams();	
	}

	/**
	 * Set the parameters for the simulation
	 * From Drouet and Pagès 2003
	 * @category architecture
	 */
	public void setParams(){
		
		plasto = (arp + (brp*getDeltaST()) + (crp*Math.pow(getDeltaST(), 2)) + (drp*Math.pow(getDeltaST(), 3))) * parms.leafApparitionRate; 
//		System.out.println(plasto);
//		System.out.println(parms.leafApparitionRate);
//		System.out.println("----------");
		// Set up the maximal number of adventive roots
		if (nodeLevel < 6) maxAdvRoot = 3;
		else if (nodeLevel==6) maxAdvRoot=4;
		else if (nodeLevel==7) maxAdvRoot=9;
		else if (nodeLevel==8) maxAdvRoot=12;	
	}

	/**
	 * Developement of the meristem, producing new leaf, stem and root
	 * @category architecture
	 */
	public void develop() {
		setParams();
		
		if(Util.isDayTime() && leafNum < 18){  
			nbWaitingPlasto += plasto;
			
			while (nbWaitingPlasto > 1.0f || firstLeaf) {
		    	
				// create and add the newly created network to the existing one
				Network stem = createNewStem();
				if(stem != null) network.addProximalNetwork(stem);
	
				relativePosition = 1.0; 
				angPhyllo += Math.PI/2;			// Phyllotaxic angle for the leaf
				angPhylloAdvRoot += Math.PI;	// Phyllotaxic angle for the roots
				spinAngle  = Math.PI;			// Spin angle
				insertionAngle = 0;				// Insertion angle for the leaf and stem
				radialAngle = 0.0;				// Radial angle for the leaf and stem
				nodeLevel = nodeLevel + 1;		// Level of the current node
				firstLeaf = false;				// Flag for the first leaf
				nbWaitingPlasto-=1; 			// Decrement of 1 
			} 
			updateThermalAge();					// Update the age
		}
	}
  
	/**
	 * Create a new stem segment, with, if required, new leaf and adventious root meristems
	 * This creation can be function of the available carbon. Priority rules for the C allocation are: 
	 * leaf meristem, stem segment, adventious root meristem
	 * @category architecture
	 * @return the newly created stem network
	 */
	public Network createNewStem() {
		// Set the available carbon for every article type
		if(parms.resolveCarbon){
			// The leaf
			leafDemand = Util.getSphereVolume(0.01) * (carbonGrowthEfficiency/2) * dryMassRatio;
			if(leafDemand < cReserve){
				leafEff = 1;
				cReserve -= leafDemand;
			}
			else{
				leafEff = cReserve / leafDemand;
				cReserve = 0;
			}
			if(cReserve < 0) cReserve = 0;
			// The stem
			stemDemand = Util.getCylinderVolume(diameter, 0.01) * (carbonGrowthEfficiency/2) * dryMassRatio;
			if(stemDemand < cReserve){
				stemEff = 1;
				cReserve -= stemDemand;
			}
			else{
				stemEff = cReserve / stemDemand;
				cReserve = 0;
			}
			if(cReserve < 0) cReserve = 0;
		}
		// If there enough C for the leaf, create the network
		if(leafEff == 1){
			
			// Create the stem segment
			Article en = new SegStem(relativePosition, insertionAngle, radialAngle, spinAngle, diameter, 
					0.01*stemEff , nodeLevel, ordre, type, parms, plasto);
			Network resEntnd=new Network(en);
	
			// Create the leaf meristem
			Article f = new MerisLeaf(0.99, insertionAngle, angPhyllo, 0, 0.01f, ordre+1, type+1, parms, leafNum, (SegStem) en);
			Network resFeuille=new Network(f);
			leafNum += 1;
			resEntnd.addDistalChildNetwork(resFeuille);
			
			// Create the adventive roots
			// From Pagès et al 1989
			if ( nodeLevel > 1 && nodeLevel <= 9 && parms.adventives) {
	    	
				for (int i=1; i<=maxAdvRoot; i=i+1 ){
					double ang = angPhylloAdvRoot-((i * Math.PI * 2.0 / maxAdvRoot) + Math.PI);
				
					// Computed the available carbon
					if(parms.resolveCarbon){
						rootDemand = Util.getSphereVolume(parms.crownRootInitDiam) * carbonGrowthEfficiency * dryMassRatio;
						if(rootDemand < cReserve){
							rootEff = 1;
							cReserve -= rootDemand;}
						else{
							rootEff = cReserve / rootDemand;
							cReserve = 0;}
						if(cReserve < 0) cReserve = 0;		
					}
					// If the root meristem is big enough to be iniated
					if(parms.crownRootInitDiam * rootEff > dr0){
						Article s =new MerisRoot(0.99,MerisRoot.angRamifAdv,ang,0.0,parms.crownRootInitDiam * rootEff, ordre+1, nodeLevel, 30, parms, Util.getNextRootID()); //0.01,angRamif,angLat,0.0,
						Network resAdvRoot = new Network(s);
						resEntnd.addDistalChildNetwork(resAdvRoot);
					}
				}
				angPhylloAdvRoot = angPhylloAdvRoot - Math.PI/6;
			}   
			return resEntnd;
		}
		else return null;
	}
 	
	/**
	 * Get the ABA consumption in this article. Based on ABA halflife
	 * From Ren et al 2007
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the aba consumption
	 */
	public double getABAConsumption() {
		double deg = 0;
		if(parms.ABASolveMethod){
			double qABA = this.envEndo.getQuantSolute(1);
			double timeStep = Time.getTimeStep();
			double halflife = ABAHalfLife; 
			deg = qABA * (1 - Math.pow(2, -timeStep / halflife));
		}
		return deg;
	}

	/**
	 * Get the quantitative production of ABA in this case, 0
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the aba production
	 */
	public double getABAProduction() {
		return 0; 
	}
	  
	/**
	 * Set the growth efficiency coefficient and the quantity of C available
	 * @category carbon
	 */
	public void setGrowthEff(double g){
		growthEfficiency = g;
		cReserve = growthDemand * g;	
	}
	
	/**
	 * Get the growth demand for this meristem, taking into account the potential formation of new
	 * stem segment and leaf and root meristems
	 * @category carbon
	 * @return the growth demand in carbon
	 */
	public double getGrowthDemand(){
		
		setParams();
		double plInit = nbWaitingPlasto;
		double demand = 0;
		boolean fLeafInit = firstLeaf;
		
		setParams();  
		if(Util.isDayTime()){  
			if(leafNum<18){	
				nbWaitingPlasto+=plasto;
			    while (nbWaitingPlasto > 1.0f || firstLeaf) { 
			    	// Stem segment
			    	demand += Util.getCylinderVolume(diameter, 0.01) * (carbonGrowthEfficiency/2) * dryMassRatio;		
			    	// Leaf meristem
			        if (ordre==1) demand += Util.getSphereVolume(0.01f) * (carbonGrowthEfficiency/2) * dryMassRatio;
			        // Root meristems
			        if ( nodeLevel > 1 && nodeLevel<=9 && parms.adventives) {
			        	for (int i = 1; i <= maxAdvRoot; i = i+1 ){
			        		demand += Util.getSphereVolume(parms.crownRootInitDiam) * carbonGrowthEfficiency * dryMassRatio;
			        	}
			    	}
			        nbWaitingPlasto-=1;
			        firstLeaf = false;
			    } 
			}
		}
		// Re Initialize the parameters
		nbWaitingPlasto = plInit;
		firstLeaf = fLeafInit;

		growthDemand = demand;
		return demand;
	}
	
	/**
	 * Get the maintenance demand for the article at this time step
	 * Used to compute the over plant requirement for maintenance
	 * From Drouet and Pagès 2003
	 * @category carbon
	 * @return the maintenance demand
	 */
	public double getMaintenanceDemand() {
		float temp = (float) (ExogenousEnvironment.getTemperature(node, stiTemp) - tRef) / 10;		
		return maintenanceRespirationRate * getVolume() * dryMassRatio * Math.pow(tempCoefficient, temp) * getDeltaST();
	}
	
	/**
	 * Get the water radial conductance of the meristem.
	 * This conductance is low, as the meristem does not lose water
	 * @category water
	 * @return the water radial conductance
	 */	
	public double getWaterRadialConductance() { return Util.epsilonDouble;}

	/**
	 * Get the water axial conductance of the meristem.
	 * This conductance is low, as the meristem as no xylem vessel yet.
	 * The meristem is therefore isolated from the rest of the root.
	 * @category water
	 * @return the water axial conductance
	 */	
	public double getWaterAxialConductance() { return Util.epsilonDouble;}

	/**
	 * Get the effect of AQP on the radial water conductivity, in this case, 1
	 * @category water
	 * @return 1
	 */
	public double getAQPCondEffect() { return 1; }

	/**
	 * Return the equivalent conductance of the article, in this case, zero
	 * @category architecture
	 * @return 0
	 */
	public double getEquivalentConductance() { return 0; }

	/**
	 * Get the article identifier of the article, in this case, the nodeLevel
	 * @category architecture
	 * @return  the nodeLevel
	 */
	public int getArticleID(){ return nodeLevel; }
	
	/**
	 * Get the growth of the article, in this case, 0
	 * @category architecture
	 * @return 0 
	 */	
	public float getGrowth() {return 0;}
	
	/**
	 * Get the distance from the base, in this case, 0
	 * @category architecture
	 * @return 0 
	 */
	public float getDistFromBase() { return (float) 0;}
	
	/**
	 * Get the percentage of radial conductivity loss, in this case, 1
	 * @category water
	 * @return  1
	 */
	public float getLrPercent() { return (float) 1; }
	
	/**
	 * Get the photosynthesis of the article, in the case of a root, 0
	 * @category carbon
	 * @return  0
	 */
	public double getPhotosynthesis() { return 0; }

	/**
	 * Get the surface of the article, in this case, a sphere
	 * @category architecture
	 * @return the surface
	 */
	public double getSurface(){ return Util.getSphereSurface(diameter); }
	
	/**
	 * Get the length of the article, in this case, 
	 * as the mersitem is a sphere, return the diameter
	 * @category architecture
	 * @return the length
	 */	  	
	public double getLength() { return diameter; } 

	/**
	 * Get the width of the article, in this case, 
	 * as the mersitem is a sphere, return the diameter
	 * @category architecture
	 * @return the width
	 */		  
	public double getWidth() { return diameter; } 

	/**
	 * Get the thickness of the article, in this case, 
	 * as the mersitem is a sphere, return the diameter
	 * @category architecture
	 * @return the thickness
	 */		  
	public double getThickness() { return diameter; } 

	/**
	 * Get the volume of the article, in this case, 
	 * as the mersitem is a sphere, volume is based on the diameter
	 * @category architecture
	 * @return the volume
	 */	
	public double getVolume() { return Util.getSphereVolume(diameter); }
	
}
