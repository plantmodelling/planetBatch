package maize;

import java.util.Random;

/**
 * This class represent a root meristem, responsible for the creation new root
 * segments and root meristem (hence the construction of the entire root system)
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium) *
 */

public class MerisRoot extends Article {
 
	static final Random random = Util.random;

	// Architecture
	public double angLat = 2.0;    						// Emission angle of the lateral roots
	public double insertAngle; 							// Insertion angle of the laterals
	public double angRamif;								// Emission angle of the ramifications
	public static double angRamifAdv = Math.PI/2;    	// Emission angle of the adventive ramifications
	     
	public double plastoRate;							// plastochronal speed used for the creation of new segments
	public double gravitropism;  						// gravitropism value of the root	
	public double growthEndAge;							// Growth end age, in thermal time									
	public double growthStartDate; 						// Age of meristem before the first lateral root initation [C.h]
	public double nbWaitingPlasto = 0.0;				// Number of plastochon waiting

	public int rootID;									// Identifier of the root (one per root)
	public int ramifLevel;								// Level of ramification
	public double diameter;       						// Diameter of the article. The meristem is condered as a sphere
	public double distanceFromBase;						// Distance from the base of the plant [cm]
	
	// Growth rate parameter based on Drouet and Pagès 2003. Growth rate is function of the root diameter
	float growthRate = 0;
	private static final float remax = 0.5f; 			// Maximal growth rate [cm/°Cd]. Init value = 0.6
	private static final int bed = 4;					// Minitial slope of the curve [°Cd]. Init value = 7
	private static final float dr0 = 0f;//.14f;			// Minimal root diameter for growth [cm]

	// Carbon
	public float dryMassRatio = 0.05f;					// Dry mass ration of the article [gC/m3]. From Drouet and Pagès 2003		
	
	// ABA
	public float ABAHalfLife = 0.7f;					// Half life of ABA [h] (Jia et al. 1996)
	
	// Positionning of the root meristem
	double count = 0.0;
	boolean counter = true;
	boolean bornInSoil = true;

	
	/**
	 * Constructor
	 */
	public MerisRoot() {
		super();
	}

	/**
	 * Constructor
	 * @param posRel
	 * @param angIns
	 * @param angRad
	 * @param angSpin
	 * @param diametre
	 * @param ordre
	 * @param type
	 * @param parms
	 * @param rid: root ID
	 */
	public MerisRoot(double posRel, double angIns, double angRad, double angSpin, 
			double diametre, int ordre, int type, ParameterSet parms, int rid) {
		super(posRel, angIns, angRad, angSpin, ordre, type, parms);
		this.diameter=diametre;
		setParams(type);
		this.rootID = rid;
	} 

	/**
	 * Constructor
	 * @param posRel
	 * @param angIns
	 * @param angRad
	 * @param angSpin
	 * @param diametre
	 * @param ordre
	 * @param node
	 * @param type
	 * @param parms
	 * @param rid
	 */
	public MerisRoot(double posRel, double angIns, double angRad, double angSpin, 
			double diametre, int ordre, int node, int type, ParameterSet parms, int rid) {
		super(posRel, angIns, angRad, angSpin, ordre, type, parms);
	    this.diameter=diametre;	
	    setParams(type, node);
	    this.rootID = rid;
	} 
  
	/**
	 * Set the parameter of the meristem based on the root type and the node (for the crown root)
	 * @category architecture
	 * @param type; the type of root
	 * @param node: the bearing node (for crown roots)
	 */
	public void setParams(int type){
		setParams(type, 0);
	}
  
	/**
	 * Set the parameter of the meristem based on the root type and the node (for the crown root)
	 * @category architecture
	 * @param type: the type of root
	 * @param node: the bearing node (for crown roots)
	 */
	public void setParams(int type, int node){
	  
		growthRate = (remax * (float) (1-Math.exp((-bed * (diameter - dr0))/remax))) / 24; // [°Cd] Drouet & Pagès 2003 
		ramifLevel = parms.branchingLevel;
	  
		if(type == 10 || type == 20 || type == 30){
			// for the need of the simulation for the PRTHH4 paper
//			growthRate = 0.01f;
			//////////////////////////////////////////////////////
			growthEndAge = parms.primaryMaxLength / growthRate;
			plastoRate = parms.primaryInterBranch / growthRate;
			growthStartDate = 0;
			angRamif = parms.secondaryInsertAngle;
		}
		
		if(type == 20){
			growthStartDate = 24;
		}
	  
		if(type == 11 || type == 21 || type == 31){
			// For the need of the simulation for the PRTHH4 paper
//			growthRate = 0.002f;
			//////////////////////////////////////////////////////
			plastoRate = parms.secondaryInterBranch / growthRate;
			growthStartDate = parms.primaryLAUZ / growthRate;
			growthEndAge =(parms.secondaryMaxLength / growthRate) + growthStartDate ;
			angRamif = parms.tertiaryInsertAngle;
			gravitropism = parms.secondaryInterBranch / parms.secondaryRootGravitropism;
		}
	  
		if(type == 12 || type == 22 || type == 32){
			growthEndAge = parms.tertiaryMaxLength / growthRate;
			plastoRate = parms.secondaryInterBranch / growthRate;
			gravitropism = parms.secondaryInterBranch / parms.tertiaryRootGravitropism;
			growthStartDate = parms.secondaryLAUZ / growthRate;
		}		
		if(type == 10) gravitropism =  parms.primaryInterBranch / parms.principalRootGravitropism;
		if(type == 20) gravitropism =  parms.primaryInterBranch / parms.seminalRootGravitropism;
		if(type == 30) gravitropism =  parms.primaryInterBranch / (4 * Math.exp(0.27 * node));
	}
  

	/**
	 * Ensures the effective development of the meristem by changing its attributes: 
	 * size, age, nbPlastoAttente, azimuthal angle and direction and 
	 * possibly generates new pieces in the network
	 * @category architecture
	 */
	public void develop() {
		  
		Network porteur=network.parentNetwork;
		
		// Define if the meristem is in the ground. If it stays too much time out of the ground it will stop growing
		if(node.z > 0) count += getDeltaST();		  
		//if (count > 200) counter = false; 

		// if the root was initiate in the soil and get out of the soil, it will stop growing
		if(bornInSoil && node.z > 0) return;
		
		// Increment the number of waiting plastochron if the segment is older than growthStartDate
		if(thermalAge > growthStartDate) nbWaitingPlasto += getDeltaST() / plastoRate;    

		// Create seminal roots (PRTHH4)
		if(Time.getTime() == 1){
//			network.addProximalNetwork(createNewSeminalRoot());	
		}
		
		
	    while (nbWaitingPlasto > 1.0f && thermalAge < growthEndAge && thermalAge > growthStartDate && counter) {// 

	    	// Create the new root
	    	network.addProximalNetwork(createNewRoot());
	    	relativePosition = 1.0;
	      
	    	insertionAngle = Util.genAleat.nextDouble() * insertAngle;		// random insertion angle
	    	radialAngle = Util.genAleat.nextDouble() * 2.0 * Math.PI; 		// random radial angle	      
	    	angLat = Util.genAleat.nextDouble() * 2.0 * Math.PI; 			// random lateral angle
	    	spinAngle = 0.0;
	      
	    	// put the newly created article in the right position
	    	porteur.position();
	      
	    	// Apply deviation to the root, mainly gravitropism
	    	float zdev = (float) (Math.max(0.01, gravitropism + gravitropism/5 * random.nextGaussian()));
	    	float xdev = (float) (parms.rootTortuosity * random.nextGaussian());
	    	float ydev = (float) (parms.rootTortuosity * random.nextGaussian());
	    	if(Util.isOutOfSoilDomain(node, soil)) zdev = 1;
	    	applyDeviation(xdev, ydev, -zdev);
	      
	    	nbWaitingPlasto -= 1; 
	    } 
	    updateThermalAge();
	} 

	/**
	 * Create a new root to insert on the network. Also create ramification (new meristem)
	 * @category architecture
	 * @return the newly created root network
	 */
	public Network createNewRoot() {

		Article s;
		Network resMetaRac = null;
		 
		// Create the new root segment
		s = new SegRoot(relativePosition, insertionAngle, radialAngle, spinAngle, diameter, ordre, type, parms, rootID);
		resMetaRac=new Network(s);
		 
		int d = type / 10;
		int n = type - (d * 10);

		// Childiameter is a function of the parent diameter
		double childDiameter = diameter/3;
		if(parms.diamStochasticity) childDiameter = (diameter / 3) + (diameter / 10) * random.nextGaussian();
		childDiameter = childDiameter * growthEfficiency;
		 
		// Create a ramification
		if (n < ramifLevel && node.z < 0 && childDiameter > dr0 ){
			s = new MerisRoot(0.01, angRamif, angLat, 0.0, childDiameter, ordre+1, type+1, parms, Util.getNextRootID());
			Network resRacSec=new Network(s);	       
			resMetaRac.addDistalChildNetwork(resRacSec);
		}
		// Create a ramification
//		if (n < ramifLevel && node.z < 0 && childDiameter > dr0 && type == 10){
//			s = new MerisRoot(0.01, angRamif, angLat, 0.0, childDiameter, ordre+1, type+1, parms, Util.getNextRootID());
//			Network resRacSec=new Network(s);	       
//			resMetaRac.addDistalChildNetwork(resRacSec);
//			for(int i = 1; i < parms.numberOfSeminals; i++){
//			s = new MerisRoot(0.01, angRamif, angLat+(i * Math.PI/3), 0.0, parms.seminalRootInitDiam/3, ordre+1, type+1, parms, Util.getNextRootID());
//			resRacSec=new Network(s);	       
//			resMetaRac.addDistalChildNetwork(resRacSec);
//			}
//		}
		
	    return resMetaRac;
	} 
	
	
	/**
	 * Create a new root to insert on the network. Also create ramification (new meristem)
	 * class created for the need of the PRTHH4 chapter
	 * @category architecture
	 * @category PRTHH4
	 * @return the newly created seminal root network
	 */
	public Network createNewSeminalRoot() {

		Article s;
		Network resMetaRac = null;
		 
		// Create the new root segment
		s = new SegRoot(relativePosition, insertionAngle, radialAngle, spinAngle, diameter, ordre, type, parms, rootID);
		resMetaRac=new Network(s);
		double inc = (2 * Math.PI) / parms.numberOfSeminals;
		for(int i = 0; i < parms.numberOfSeminals; i++){
			// Create a seminal root
			double stochasticity = 1;
			if(parms.semAngStochasticity) stochasticity =  Util.genAleat.nextDouble();
			double ang1 = 1.5 * angRamif;
			double ang2 = (angLat + (i * inc)) * stochasticity;
			s = new MerisRoot(0.1, ang1, ang2, 1.0, parms.seminalRootInitDiam, ordre+1, 20, parms, Util.getNextRootID());
			Network resRacSec=new Network(s);	       
			resMetaRac.addDistalChildNetwork(resRacSec);
		}
		return resMetaRac;
	} 
	
	/**
	 * Get the growth demand for this meristem, taking into account the potential formation of new
	 * root segments and meristem
	 * @category carbon
	 * @return the growth demand in carbon
	 */
	public double getGrowthDemand(){
		
		double nbWaitingPlastoInit = nbWaitingPlasto;
		double growthDemand = 0;
		int d = type / 10;
		int n = type - (d * 10);
		
		if(thermalAge > growthStartDate) nbWaitingPlasto += getDeltaST() / plastoRate;    
	    while (nbWaitingPlasto > 1.0f && thermalAge < growthEndAge && thermalAge > growthStartDate && counter) {// 
	    	growthDemand += carbonGrowthEfficiency * Util.getCylinderVolume(diameter, diameter / 10) * dryMassRatio;
			if (n < ramifLevel && node.z < 0 ){
				growthDemand += Util.getSphereVolume(diameter / 3) * carbonGrowthEfficiency * dryMassRatio;
			}
	    	nbWaitingPlasto -= 1; 
	    } 
	    // re-initialise nbWaitingPlasto
	    nbWaitingPlasto = nbWaitingPlastoInit;
		return growthDemand;
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
			double timeStep = Time.getTimeStep() / parms.resolveIteration;
			double halflife = ABAHalfLife;
			deg = qABA*(1-Math.pow(2, -timeStep / halflife));
		}
		return deg;
	}

	/**
	 * Get the quantitative production of ABA, based on the volume and the time step
	 * From Dodd 2008
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the aba production
	 */
	public double getABAProduction() {
		double pseudoWP = getWaterPotExo();
		double prod = 0.000005 * (-170 * pseudoWP) * Time.getTimeStep() * getVolume(); 
		return prod;
	} 

	/**
	 * Set the distance between the article and the base of the root it is being part of
	 * @category architecture
	 */
	public void setDistanceFromBase(){
		SegRoot p;
		if(getParentArticle() != null){
			if(getParentArticle() instanceof SegRoot){
				p = (SegRoot) getParentArticle();
				if(p.getArticleID() == getArticleID()) {
					p.setDistanceFromBase();
					distanceFromBase = p.getDistFromBase() + getLength();
				}
				else distanceFromBase = getLength();
			}	  
		}
	}
 
	/**
	 * Get the water radial conductance of the meristem.
	 * This conductance is high, as the meristem as no secondary structure yet.
	 * The conductance is low if the meristem is in the air
	 * @category water
	 * @return the water radial conductance
	 */
	public double getWaterRadialConductance() { 
		if(node.z > 0) return Util.epsilonDouble;
		else return 1 / Util.epsilonFloat;
		}

	/**
	 * Get the water axial conductance of the meristem.
	 * This conductance is low, as the meristem as no xylem vessel yet.
	 * The meristem is therefore isolated from the rest of the root.
	 * @category water
	 * @return the water axial conductance
	 */
	public double getWaterAxialConductance() { return Util.epsilonDouble; }
	  
	/**
	 * Return the equivalent conductance of the article, in this case, zero
	 * @category architecture
	 * @return 0
	 */
	public double getEquivalentConductance() { return 0; }

	/**
	 * Get the article identifier, in this case, the rootID
	 * @category architecture
	 * @return the rootID
	 */
	public int getArticleID(){ return rootID; }


	/**
	 * Get the distance from the base of the root
	 * @category architecture
	 * @return the distance from the base of the root
	 */
	public float getDistFromBase(){return (float) distanceFromBase;}

	/**
	 * Get the growth of the article. Zero for the meristem
	 * @category architecture
	 * @return 0	 
	 s*/
	public float getGrowth() { return 0; }

	/**
	 * Get the effect of AQP on the radial water conductivity, in this case, 1
	 * @category water
	 * @return 1
	 */
	public double getAQPCondEffect() { return 1; }

	/**
	 * Get the reduction in radial water conductivity, in this case, 1
	 * @category water
	 * @return 1
	 */
	public float getLrPercent() { return 0; }

	/**
	 * Geth the photosynthesis of the article, in this case, 0
	 * @category carbon
	 * @return 0
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
