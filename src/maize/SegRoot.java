package maize;

import java.lang.Math;

/**
 * This class represent a root segment.
 * Its aims is to build the structure of the plant.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium) 
 */

public class SegRoot extends Article {

    // Developement
	public double growthStartAge = 0.0;  			// Growth start age, in thermal time
	public double growthEndAge;  					// Growth end age, in thermal time
	
	// Architecture
	public int rootID;								// Root identifier
	public double maxLength;						// Maximal length of the root segment
	public double distanceFromBase;					// Distance if this root segment to the base of the root
  
	// ABA
	public float ABAHalfLife = 0.7f;				// Half life of ABA [h] (Jia et al. 1996)
	
	// Growth rate parameter based on Drouet and Pagès 2003. Growth rate is function of the root diameter
	public float axialGrowthRate; 
	private static final float remax = 0.5f; 			// Maximal growth rate [cm/°Cd]. Init value = 0.6
	private static final int bed = 4;					// Minimal slope of the curve [°Cd]. Init value = 7
	private static final float dr0 = 0.0f;//14f;		// Minimal root diameter for growth [cm]

	// Parameters for the final dry mass of the stem segment(from Drouet 2003)
	private float rmvMax = 0.14f;					// Maximal dry mass ratio [gC/cm3]									
	private float rmvMin = 0.05f;					// Minimal dry mass ratio [gC/cm3]
	private float dryMassRatio = rmvMin;			// Actual dry mass ratio [gC/cm3]
	private double dryMassRatioGrowth = 0;			// Growth of the dry mass ratio

 
	/**
	 * Constructor
	 */
	public SegRoot(){
		super();
	}
  
	/**
	 * Constructor
	 * @param posRel
	 * @param angIns
	 * @param angRad
	 * @param angSpin
	 * @param diamInit
	 * @param ordre
	 * @param t: type
	 * @param parms
	 * @param rid
	 */
	public SegRoot(double posRel, double angIns, double angRad, double angSpin, 
			double diamInit, int ordre, int t, ParameterSet parms, int rid) {
		super(posRel,angIns,angRad,angSpin, ordre, t, parms); 
		diameter=diamInit;  
		length=diamInit / 10;
		setParams(type);
		rootID = rid;
		type = t;
	}
  
	/**
	 * Set the root segment parameter for growth based on the root segment type
	 * @category architecture
	 */
	public void setParams(int type){
	  
		axialGrowthRate = (remax * (float) (1-Math.exp((-bed*(this.diameter-dr0))/remax))) / 24; // [°Cd] Drouet & Pagès 2003 

		if(type == 10 || type == 20 || type == 30){
			// for the need of the simulation for the PRTHH4 paper
//			axialGrowthRate = 0.01f;
			//////////////////////////////////////////////////////
			maxLength = parms.primaryInterBranch;
		}
		else {
			// for the need of the simulation for the PRTHH4 paper
//			axialGrowthRate = 0.002f;
			//////////////////////////////////////////////////////
			maxLength = parms.secondaryInterBranch;
		}
	
		growthEndAge = maxLength / axialGrowthRate;  
		dryMassRatioGrowth = (rmvMax - rmvMin) / growthEndAge;
	}
  
	/**
	 * Insure the growth and developement of the root segment
	 * @category architecture
	 */
	public void develop() {

		if(thermalAge > growthStartAge && thermalAge < 1.5 * growthEndAge){

			// Increase the length
			length += getAxialGrowth() * growthEfficiency;
			if(length > maxLength) length = maxLength;
		
			// Increase the dry mass
			if(dryMassRatio < rmvMax) dryMassRatio += dryMassRatioGrowth;
			if(dryMassRatio > rmvMax) dryMassRatio = rmvMax;
			dryMass = getVolume() * dryMassRatio * growthEfficiency;	        
		}
	    updateThermalAge(); 
	    
	    // Does not allow the segment to have a water potential value below the permanent wilting point
	    //if(getWaterPot() < parms.permanentWiltingPoint) killArticle();
	}
  
	/**
	 * Compute the growth in length of the article based on
	 * its endogenous and exogenous environment
	 * @category architecture 
	 * @return the growth in length
	 */
	public double getAxialGrowth() {

		if(thermalAge > growthStartAge && thermalAge < 1.5 * growthEndAge){  	    	   	
			return axialGrowthRate * getDeltaST() ;// * getSoilMechanicalGrowthEffect() * getABAGrowthEffect();
		}
		else return 0.0; 
	}

	/**
	 * Get the mechanical of the soil on the root growth
	 * @category architecture	
	 * @return the mechanical effect of the soil [0-1]
	 */
	public double getSoilMechanicalGrowthEffect(){
		return ExogenousEnvironment.getMechanicalConstrain(node, stiRM.setCoordinates(this.node.x, this.node.y, this.node.z));
	}
	
	/**
	 * Get the effect of ABA concentration on the growth of the root segment
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the effect of ABA on growth  [0-1]
	 */
	public double getABAGrowthEffect(){
		double aba = getConcSolute(1);
		double effetABA = 1;
		// growth of primary roots is stimulated
		if (ordre == 1){ 						
			effetABA = (500 + 3 * aba) / (500 + aba); 
		}
		// Growth of secondary and tertiary root is stopped (simplification)
		else if (ordre > 1 && aba > 150){ 			
			effetABA = 0;
		}
		return effetABA;
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
		double prod = 0;
		if(parms.ABAProduction)  prod =  9.5635 * Math.pow(Math.E, (-5.224 * getWaterPot())) * getVolume();
		return prod;
	}
	
	/**
	 * Calculate and return the axial conductance to water flux of the segment based on its age
	 * Kh is the axial root conductivity [m4/s.MPa]
	 * From Doussan 1998
	 * @category water
	 * @return axial conductance
	 */
	public double getWaterAxialConductance() {
		  
		double age = getAge() / 24;	  
		double cond = 0.0; 
		  
		// If the root is an axe
		if(type == 10 || type == 20 || type == 30){
			  
			float d1 = 12 / axialGrowthRate; 
			float d2 = 35 / axialGrowthRate;
			float d3 = 50 / axialGrowthRate;
			  		  
			if(thermalAge<=d1) cond  = Util.epsilonDouble; // 0
 			else if (thermalAge>d1 && thermalAge<=d2) cond  = (0.0147 * (thermalAge * axialGrowthRate) - 0.0146) * 1e-9;
			else if (thermalAge>d2 && thermalAge<=d3) cond  = (0.3 * (thermalAge * axialGrowthRate) - 10) * 1e-9;
			else cond =5e-9;
//			cond  = 1e-9;
		}
		  
		// If the root is a branch
		else {
			float d1 = 10; 
			float d2 = 13;
			float d3 = 20;
			float d4 = 25;
			  
			if(age<=d1) cond  = Util.epsilonDouble; // 0
			else if (age>d1 && age<=d2) cond  = (0.267 * age - 2.67) * 1e-12; //(0.0667 * age - 0.0666) * 1e-12; 
			else if (age>d2 && age<=d3) cond  = 0.8e-12;
			else if (age>d3 && age<=d4) cond  = (0.24 * age - 4) * 1e-12;
			else cond = 2e-12;
			//cond  = 1e-12;
		}	  
		 		
		return cond * getCavitationEffect() * parms.axMod;
	}

	/**
	 * Calculate and return the radial conductance to water flux of the segment based on its age
	 * Lr is the radial root conductivity [m/s.MPa]
	 * From Doussan 1998
	 * @category water
	 * @return radial conductance
	 */
	public double getWaterRadialConductance(){
		  	  
		double cond = 0.0; 
		
		// If the root is an axe
		if(type == 10 || type == 20 || type == 30){
			  
			float d1 = 12 / axialGrowthRate;
			float d2 = 20 / axialGrowthRate;
			float d3 = 45 / axialGrowthRate;
			float d4 = 60 / axialGrowthRate;
			  
			if(thermalAge <= d1) cond = 2.0e-7;
			else if (thermalAge>d1 && thermalAge<=d2) cond = (-0.17 * (thermalAge * axialGrowthRate) + 4.3) * 1e-7;
			else if (thermalAge>d2 && thermalAge<=d3) cond = 0.8e-7;
			else if (thermalAge>d3 && thermalAge<=d4) cond = (-0.04 * (thermalAge * axialGrowthRate) + 2.6)* 1e-7;
			else cond=0.2e-7;
//			cond=0.2e-7;
		}
		  
		// If the root is a branch
		else {
			double age = getAge() / 24;	  

			float d1 = 10;
			float d2 = 16;
			  
			if(age <= d1) cond = 2.0e-7;
			else if (age>d1 && age<=d2) cond = (-0.3 * age + 5) * 1e-7;
			else cond=0.2e-7;
//			cond=0.2e-7;
		}
		   		  
		// If the segment is in the air, its resistance is too big to uptake or lose water
		if(node.z > 0) cond = Util.epsilonDouble;
		  
		return cond * getAQPCondEffect() * getABACondEffect() * parms.radMod;
	}

	/**
	 * Logisitic function between 1 and 10 reflecting the influence of the aquaporins as a function of the water potential
	 * If WP < 0.6, mod = 1 and if WP > 1.2 = 10;
	 * @category water
	 * @return the AQP effect
	 */
	public double getAQPCondEffect(){
		float A = 0.4f;
		int K = 1;
		float B = 8f;
		float Q = 0.3f;
		float v = 0.3f;
		float m = -1.1f;
		double x = getWaterPot() + parms.aqpMod;
		if(x > 0) x = 0;
		
		if(parms.aqp) return A + ((K - A) / Math.pow((1 + (Q * Math.exp(-B * (x-m)))),(1/v)));
		else return 1;
	}
	
	/**
	 * The modificator of ABA is calculated on the basis of these data coming from Hose(2000) : 10nM -> 1x; 100nM -> 3x; 1000nM -> 4x.
	 * It is calculated by a logistic function : ((A-D)/1+(Math.pow(Math.log(getConcSolute(1))/C,B)))+D
	 * Where : A = minimum; B = Slope factor; C = inflexion point; D = maximum
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the effect of ABA on the root radial conductivity
	 */
	public double getABACondEffect(){
		
		if(getConcSolute(1)  >= 1){
			return ((1 - 4) / 1 + (Math.pow(Math.log(getConcSolute(1)) / 4, 5.5))) + 4;
		}
		else return 1;		   		
	}
	
	/**
	 * Get the demand in carbon for the growth of the article.
	 * Takes into account the growth in length of the article
	 * From Drouet and Pagès 2003
	 * @category carbon
	 * @return the growth demand in C
	 */
	public double getGrowthDemand(){
		return carbonGrowthEfficiency * Util.getCylinderVolume(diameter, getAxialGrowth()) * dryMassRatio; 
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
		return (1-parms.aerenchyma) * maintenanceRespirationRate * getVolume() * dryMassRatio * Math.pow(tempCoefficient, temp) * getDeltaST();
	}
	
	/**
	 * Compute the distance between the article and the base of its bearing root
	 * @category architecture
	 */
	public void setDistanceFromBase(){
		SegRoot p;
		if(getParentArticle() != null){
			if(getParentArticle() instanceof SegRoot){
				p = (SegRoot) getParentArticle();
				if(p.getArticleID() == getArticleID()){
					p.setDistanceFromBase();
					distanceFromBase = p.getDistFromBase() + getLength();
				}
				else distanceFromBase = getLength();	
			}	  
		}
	}

  
    /**
     * Get the equivalent conductance of the segment
     * @author Guillaume Pagès
     * @category water
     * @return the equivalent resistance
     */
	public double getEquivalentConductance() {
		double s = getSurface() / 1.0e4 * 1 / getWaterRadialResistance();
		if (network.childNetworkList != null) {
			for (Network resFils : network.childNetworkList) {
				double c1 =  1 / resFils.baseArticle.getWaterAxialResistance() / resFils.baseArticle.getLength() / 1.0e2;
				double c2 = resFils.baseArticle.getEquivalentConductance();
				if ( c1 + c2 != 0) s += (c1 * c2) / (c1 + c2);
			}
		}
		KhEqu = s;
		return s;	 
	}
	
	/**
	 * Get the article identifier of the article, in this case, the rootID
	 * @category architecture
	 * @return  the rootID
	 */
	public int getArticleID(){ return rootID; }
	
	/**
	 * Get the growth of the article
	 * @category architecture
	 * @return the growth rate 
	 */	
	public float getGrowth() {return (float) getAxialGrowth();}
	
	/**
	 * Get the distance between the article the base of the root
	 * @category architecture
	 * @return the distance from the base of the root 
	 */
	public float getDistFromBase() { return (float) distanceFromBase;}
	
	/**
	 * Get the percentage of radial conductivity loss, in this case, the effect
	 * of aquaporins
	 * @category water
	 * @return  the effect of AQP
	 */
	public float getLrPercent() { return (float) getAQPCondEffect(); }
	
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
