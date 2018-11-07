package maize;

/**
 * This class represent a stem segment.
 * Its aims is to build the structure of the plant.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium) 
 */

public class SegStem extends Article {

  
	public double growthStartAge = 0.0;  	// Growth start age, in thermal time
	public double growthEndAge = 10000.0;  	// Growth end age, in thermal time
	public int nodeLevel;					// Level of the node
	
	// Water
	public float xylemPercentage = 0.017f;	// Percentage of the stem which is xylem vessels
	
	// Parameters for the final length of the stem segment(from Drouet and Pagès 2003)
	public int n0 = 5;
	public float aLi = 0.202f;
	public float bLi = 1.233f;
	public float cLi = -0.04f;
	public int LiM = 23;
	public float stemGrowthRate = 0.41f;
	public float finalLength ;				// Final length of the article, depend on the node level
	public float axialGrowthRate;			// Axial growth rate of the article 
  
	// Parameters for the final diameter of the stem segment(from Drouet and Pagès 2003)
	private double abi = 2.3, bbi = 3.3, cbi = -0.17f;
	private int nbi = 7;
	private double finalDiameter;			// Final diameter of the article, depend on the node level
	private double initialDiameter;			// Initial diameter of the article
	private double radialGrowthRate;		// Radial growth rate of the article

	// Parameters for the final dry mass of the stem segment(from Drouet and Pagès 2003)
	private float imvMax = 0.14f;//0.01f;			// Max dry mass ratio [gC/cm3]
	private float imvMin = 0.005f;			// Min dry mass ratio [gC/cm3]
	private double dryMassRatioGrowth;		// Growth of the dry mass ration
	private double dryMassRatio = imvMin;	// Current dry mass ratio [gC/cm3]	

	// ABA
	public float ABAHalfLife = 0.7f;		// Half life of ABA [h] (Ren et al. 2007)
	
	/**
	 * Constructor
	 */
    public SegStem() {
    	super();
    }

    /**
     * Constructor
     * @param posRel
     * @param angIns
     * @param angRad
     * @param angSpin
     * @param diam
     * @param l
     * @param nodeCount
     * @param ordre
     * @param type
     * @param parms
     */
    public SegStem(double posRel, double angIns, double angRad, double angSpin, 
		  double diam, double l, int nodeCount, int ordre,int type, ParameterSet parms) {
   
    	super(posRel,angIns,angRad,angSpin, ordre, type, parms); 
    	diameter=diam;  
    	initialDiameter = diam;
    	length = l;
    	nodeLevel=nodeCount;
    	setParams();
    }

    /**
     * Constructor
     * @param posRel
     * @param angIns
     * @param angRad
     * @param angSpin
     * @param diam
     * @param l
     * @param nodeCount
     * @param ordre
     * @param type
     * @param parms
     * @param plasto
     */
    public SegStem(double posRel, double angIns, double angRad, double angSpin, 
		  double diam, double l, int nodeCount, int ordre,int type, ParameterSet parms, double plasto) {
   
    	super(posRel,angIns,angRad,angSpin, ordre, type, parms); 
    	diameter=diam;  
    	initialDiameter = diam;
    	length = l;
    	nodeLevel=nodeCount;
    	setParams();
    	growthStartAge = plasto;
    }	
  
  /**
   * Set the parameters for the growth and final length of the stem segment
   * from Drouet 2003
   * @category architecture
   */
    public void setParams(){
	  
    	// Growth rate
    	axialGrowthRate = stemGrowthRate * parms.leafGrowthRate;
	  
    	// Final length of the internode which depend on the node number of the segment
    	if(nodeLevel <= n0) finalLength = 0;
    	else if(nodeLevel > n0 && nodeLevel <= n0 + 5) finalLength = LiM * (aLi * (nodeLevel - n0));
    	else finalLength = LiM * (bLi + (cLi * (nodeLevel - n0)));
    	if(finalLength <= 0) finalLength = 0.01f; 
	  
    	// Growth duration
    	growthEndAge = finalLength / axialGrowthRate;
	   
    	// Final diameter
    	if (nodeLevel < nbi) finalDiameter = abi;
    	else finalDiameter = bbi+cbi * nodeLevel;
	  
    	// Dry mass and radial growth
    	if(nodeLevel > n0){
    		dryMassRatioGrowth = (imvMax - imvMin) / growthEndAge;
    		radialGrowthRate = (finalDiameter - initialDiameter) / (growthEndAge * 10) ;
    	}
    	else{
    		radialGrowthRate = 0;
    		dryMassRatioGrowth = 0;
    	}  
    }
 

    /**
     * Insure the axial and radial growth of the stem segment
     * Insure the dry mass growth
     * @category architecture
     */
    public void develop() {
    	if(Util.isDayTime()){  	
  		
    		// Axial growth
    		length +=  getAxialGrowth() * growthEfficiency;
    		if(length > finalLength) length = finalLength;
  	    
    		// Radial growth
    		diameter += getRadialGrowth() * growthEfficiency;
    		if(diameter > finalDiameter) diameter = finalDiameter;
  	    
    		// Dry mass growth
    		if(dryMassRatio < imvMax) dryMassRatio += dryMassRatioGrowth;
    		if(dryMassRatio > imvMax) dryMassRatio = imvMax;	    
    		dryMass = getVolume() * dryMassRatio * growthEfficiency;
  	    
    		// Age growth
    		updateThermalAge();
    	}
    }

    /**
     * Get the axial growth of the article
     * @category architecture
     * @return the axial growth
     */
    public double getAxialGrowth() {
    	if ((thermalAge >= growthStartAge) && (length < finalLength)){
    		return axialGrowthRate * getDeltaST();
    	}
    	else return 0.0f;
    }

    /**
     * Get the radial growth of the article
     * @category architecture
     * @return the radial growth
     */
    public double getRadialGrowth() {

    	if ((thermalAge >= growthStartAge) && (diameter < finalDiameter)) { 
    		return radialGrowthRate * getDeltaST();
    	}
    	else return 0.0;
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
   * Axial conductance is based on Poisseuille Law
   * @category water
   * @see getWaterAxialResistance()
   * @return the water axial conductance
   */
    public double getWaterAxialConductance(){ return 1 / getWaterAxialResistance(); }
    
    /**
     * Axial resistance of the stem
     * Based on Poiseuille'law
     * From Li et al 2009
     * @category water
     * @return the water axial resistance
     */
    public double getWaterAxialResistance(){
    	
    	double xylemDiametre = diameter * xylemPercentage ;
    	float waterViscosity = 1.002f;
    	double axialRes = (8 * waterViscosity) / (Math.PI * Math.pow(xylemDiametre, 4));
    	    	    	
    	return axialRes * (1 / getCavitationEffect());
    }
  
    /**
     * Get the water radial conductance of the stem, which is low
     * @category water
     * @return the radial conductance
     */
    public double getWaterRadialConductance(){ return Util.epsilonDouble; }
   

    /**
     * Get the equivalent conductance of the segment
     * @author Guillaume Pagès
     * @category water
     * @return the equivalent resistance
     */
    public double getEquivalentConductance() {
    	double s = (getSurface() / 1.0e4) * (1 / getWaterRadialResistance());
    	if (network.childNetworkList != null) {
    		for (Network resFils : network.childNetworkList) {
				double c1=  1 / resFils.baseArticle.getWaterAxialResistance() / (resFils.baseArticle.getLength() / 1.0e2);
				double c2= resFils.baseArticle.getEquivalentConductance();
				if ( (c1 + c2) != 0) s += (c1 * c2) / (c1 + c2);				
			}
		}
		KhEqu = s;
		return s;	
	}

	/**
	 * Get the growth demand of the article
	 * From Drouet and Pagès 2003
	 * @category carbon
	 * @return the growth demand [gC]
	 */	
	public double getGrowthDemand(){ 
		
		if(Util.isDayTime()) return (carbonGrowthEfficiency/2) * Util.getCylinderVolume(getRadialGrowth(), getAxialGrowth()) * dryMassRatio;
		else return 0;
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
	 * Get the effect of AQP on the radial water conductivity, in this case, 1
	 * @category water
	 * @return 1
	 */
	public double getAQPCondEffect() { return 1; }

	/**
	 * Get the article identifier of the article, in this case, the rootID
	 * @category architecture
	 * @return  the rootID
	 */
	public int getArticleID(){ return nodeLevel; }
	
	/**
	 * Get the growth of the article
	 * @category architecture
	 * @return the growth rate 
	 */	
	public float getGrowth() {return axialGrowthRate;}
	
	/**
	 * Get the distance between the article the base of the root
	 * @category architecture
	 * @return the distance from the base of the root 
	 */
	public float getDistFromBase() { return (float) 0;}
	
	/**
	 * Get the percentage of radial conductivity loss, in this case 1
	 * @category water
	 * @return 1
	 */
	public float getLrPercent() { return 1; }
	
	/**
	 * Get the photosynthesis of the article, in the case of a root, 0
	 * @category carbon
	 * @return  0
	 */
	public double getPhotosynthesis() { return 0; }

	/**
	 * Get the length of the article
	 * @category architecture
	 * @return the length
	 */	  	
	public double getLength() { return length; } 

	/**
	 * Get the width of the article, in this case, 
	 * as the segment is a cylinder, return the diameter
	 * @category architecture
	 * @return the width
	 */		  
	public double getWidth() { return diameter; } 

	/**
	 * Get the width of the article, in this case, 
	 * as the segment is a cylinder, return the diameter
	 * @category architecture
	 * @return the thickness
	 */		  
	public double getThickness() { return diameter; } 

	/**
	 * Get the volume of the article, in this case, 
	 * as the mersitem is a cylinder, volume is based on the diameter and length
	 * @category architecture
	 * @return the volume
	 */	
	public double getVolume() { return Util.getCylinderVolume(diameter, length); }
	
	/**
	 * Get the surface of the article, in this case, a cylinder
	 * @category architecture
	 * @return the surface
	 */
	public double getSurface(){ return Util.getCylinderSurface(diameter, length); }
	
}
