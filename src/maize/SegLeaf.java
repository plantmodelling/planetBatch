package maize;

import java.lang.Math;

/**
 * This class represent a leaf segment.
 * Its aims is to build the structure of the plant.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium) 
 */

public class SegLeaf extends Article {
  
	// Architecture
	public double width;  								// Width of the segment	 
	public double thickness;  		 					// Thickness of the segment
	public double length;  								// Length of the segment
	public boolean lamina;								// Flag for the lamina / sheath
	private double maxLength;  							// Maximal length of the segment
	public float leafGrowthRate;  						// Growth rate of this leaf segment
	public int leafNum;									// Number of the leaf to which this semgent belong
	public double growthEndAge;							// Growth end age, in thermal time  	
	public SegStem stem;								// Stem segment bearing the leaf
	
	// Water
	public double gsMax = 1e-9;							// Maximal openning of the stomata (default = 3e-9)
  
	// ABA
	public float ABAHalfLife = 0.96f;					// Half life of ABA [h] (Ren et al. 2007)
	
	// Carbon. From Drouet and Pagès 2003
	private float dryMassRatio = 0.006f;				// Dry mass ratio [g/m3]
	private double dryMassRatioGrowth = 0;				// Growth rate of the dry mass ratio (aka tissue devellopement)
	private float lmaMin = 0.001f;//0.005f;						// Min dry mass ratio
	private float lmaMax = 0.0075f;//0.14f;						// Max dry mass ratio
	public double photosynthesis = 0;					// Photosynthesis of the leaf segmtn
	public float Vmax = parms.VMax; 
	public float Jmax = parms.JMax;
	// Carbon. From Earl 1998
	public float maintenanceRespirationRate = 7.2e-6f; 	// gC/gDW s

	/**
	 * Constructor
	 */
	public SegLeaf() {
		super(); 
	} 

	/**
	 * Constructor
	 * @param posRel
	 * @param angIns
	 * @param angRad
	 * @param angSpin
	 * @param w: width of the segment
	 * @param le: length of the segment
	 * @param t: thickness of the segment
	 * @param ordre
	 * @param type
	 * @param parms: Parameter set of the simulation
	 * @param l: leaf number of the segment
	 * @param lam: if lamina = true
	 * @param ml: maximal length	
	 * @param s: stem segment bearing the leaf
	 */
	public SegLeaf(double posRel, double angIns, double angRad, double angSpin, double w, double le, 
		  double t, int ordre, int type, ParameterSet parms, int l, boolean lam, double ml, SegStem s) {

		super(posRel,angIns,angRad,angSpin, ordre, type, parms);
	    width=w;  
	    length=le;  
	    thickness=t;      
	    lamina = lam;
	    leafNum = l;
	    maxLength = ml;
	    stem = s;
	    setParams();
	  } 
  
  
	/**
	 * Set the parameters for this leaf segment
	 */
	public void setParams(){
		leafGrowthRate = parms.leafGrowthRate;
		growthEndAge = maxLength / leafGrowthRate;
		dryMassRatioGrowth = (lmaMax - lmaMin) / growthEndAge;
	}
	
	/**
	 * Develop the leaf segment, which is, increase its length and its mass
	 * @category architecture
	 */
	public void develop() {
		if(Util.isDayTime()){
			// Increase length
			length += getAxialGrowth() * growthEfficiency;	
			if(length > maxLength) length = maxLength;
		  
			// Increase the dry mass
			if(dryMassRatio < lmaMax) dryMassRatio += dryMassRatioGrowth;
			if(dryMassRatio > lmaMax) dryMassRatio = lmaMax;
			dryMass = getVolume() * dryMassRatio * growthEfficiency;
			
			// Increase the age
			updateThermalAge();
		}    
	}

	/**
	 * Get the ABA consumption by the article
	 * ABA degradation is based on a degradation rate, documented in Ren et al. 2007 (Dynamic analysis of ABA accumulation)
	 * ABA Half-life in stress condition : 1.1066 h / non-stress condition : 0.8119 h (mean value : 0,96h)
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the ABA consumption
	 */
	public double getABAConsumption() {
		double deg = 0;
		//if(parms.ABASolveMethod){
			double qABA = this.envEndo.getQuantSolute(1);
			double timeStep = Time.getTimeStep() / parms.resolveIteration;
			deg = qABA * (1 - Math.pow(2, -timeStep / ABAHalfLife));
			//Util.log("aba = "+qABA);
			//Util.log("time = "+timeStep);
			//Util.log("deg = "+deg);
			//}
		return deg;	
	}
 
	/**
	 * Get the ABA production by the article over the time step
	 * From Dodd 2008
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the ABA production
	 */
	public double getABAProduction() {
		double prod = 0;	  
		if(parms.ABAProduction) prod = 9.5635 * Math.pow(Math.E, (-5.224 * getWaterPot())) * 0.5 * getVolume();
		//Util.log("prod = "+prod);
		return prod;
	}

  
	/**
	 * Get the growth in length of the article
	 * @category architecture
	 * @return the growth in length
	 */
	public double getAxialGrowth() {

		if (length < maxLength) {
			double effetA=1.0;
			return effetA * leafGrowthRate * getDeltaST();
		}
		else return 0.0; 
	}
  
  
	/**
	 * Get the water radial conductance of the leaf segment (aka stomata openning)
	 * This value is modified by both endogenous and exogenous parameters
	 * @category water
	 * @return the water radial conductance
	 */
	public double getWaterRadialConductance(){	  	
	  	  
	  	double gs = gsMax *  getABAModifier() * getEnviModifier() * getStomataModifier();   	
	  	if(node.z < 0) gs = 0;
//	  	if(!lamina || node.z < 0) gs = 0;
		return gs; 
	}

	/**
	 * Get stomata modifier based on leaf water potential (Cochard02) 
	 * @category water
	 * @return the stomata modifier
	 */
	public double getStomataModifier(){
	  	  
//		int A = 0;
//		int K = 1;
//		int B = 7;
//		int Q = 20;
//		float v = 0.7f;
//		float m = 0.5f;
//		double x = -getWaterPot() - parms.stomataMod;
//		if(x < 0) x = 0;
//	  			  	  
//		if(parms.stomata) return 1 - (A + ((K - A) / Math.pow((1 + (Q * Math.exp(-B * (x-m)))),(1/v))));
//		else return 1;
		
		float[] parm = new float[]{5f, 1.1f};
		
		double x = -getWaterPot() - parms.stomataMod;
		if(x < 0) x = 0;
		
//		System.out.println("stox = "+x);
		
		if(parms.stomata) return Math.exp(-Math.pow((x/parm[1]),parm[0]));
		else return 1;
		
	}
  
  
	/**
	 * Get stomata modifier based environment and on Jarvis equations (Jarvis 76) 
	 * @category water
	 * @return the stomata modifier
	 */
	public double getEnviModifier(){
	  
		int TL = 10;
	  	int TH = 50;
	  	int a1 = 502;
	  	int a2 = 28;
	  	float a3 = 0.03f;
	  	
	  	double vpdMod = Math.exp( -a3 * ExogenousEnvironment.getVPD(node)) ;
	  	double lightMod = ExogenousEnvironment.getLight(node, parms.maxPAR) / a1;
  		double tempExp = (TH-a2)/(a2-TL);
  		double temp = ExogenousEnvironment.getTemperature(node, stiTemp);
	  	double tempMod = Math.pow(((temp-TL)*(TH-temp)), tempExp) / (Math.pow(((a2-TL)*(TH-a2)), tempExp)); 
	  		
	  	return vpdMod * lightMod * tempMod;
	}
    
	/**
	 * Get stomata modifier based on leaf ABA content and water potential (Tardieu and Davies )
	 * @author Vincent Larondelle
	 * @category water
	 * @return the stomata modifier
	 */
	public double getABAModifier(){ 
	
		if(parms.ABASolveMethod){
	  		double gsMin = 0.5e-10; 			// Approximation from graphs in Tardieu & Davies. Has to be fitted
		  	double alpha = gsMax - gsMin;
		  	double beta = -4.2502e-3; 			// From Tardieu & Davies : -2.69 x 0.0157, found by recursivity
		  	double delta = 2.7770982; 			// From Tardieu & Davies : 0.183 x 15.1754, found by recursivity
		  	double aba = getConcSolute(1);
		  	double WP = -getWaterPot();//(-0.6-getWaterPot())-0.6; 
		  	// Vincent Larondelle Impossible d'expliquer pourquoi c'est comme ça. Très bizarre,
		  	// mais c'est la seule solution pour correspondre aux graphs.
		  	double mod = (gsMin + alpha * Math.exp(aba * beta * Math.exp(delta * WP))) / gsMax;
		  	if(mod * 4 > 1) return 1;
		  	else return mod * 4;
		}
		else return 1;
	}
  
	/**
	 * Get the axial conductance of the leaf segment.
	 * Based on the axial cnductance of the bearing stem segment.
	 * Modified by caviation.
	 * @category water
	 * @return the water axial conductance
	 */
    public double getWaterAxialConductance(){     	
    	return (stem.getWaterAxialConductance() / 10 ) * getCavitationEffect(); 
    }
  
    /**
     * Get the equivalent conductance of the segment
     * @author Guillaume Pagès
     * @category water
     * @return the equivalent resistance
     */
	public double getEquivalentConductance() {
		double s = (getSurface() / 1.0e4) * getWaterRadialConductance();
		if (network.childNetworkList != null) {
			for (Network resFils : network.childNetworkList) {
				double c1=  resFils.baseArticle.getWaterAxialConductance() / (resFils.baseArticle.getLength() / 1.0e2);
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
		if(Util.isDayTime()) return getAxialGrowth() * width * thickness * (carbonGrowthEfficiency/2) * dryMassRatio;	
		else return 0;
	}
	
	/**
	 * Get the maintenance demand for the article
	 * From Earl et al. 1998
	 * @category carbon
	 * @return the maintenance demand
	 */
	public double getMaintenanceDemand(){
		return maintenanceRespirationRate * dryMass * Time.getTimeStep();
	}
	
	/**
	 * Get the article identifier. for leaf, it is the leaf number
	 * @category architecture
	 * @return leaf number
	 */
	public int getArticleID(){return leafNum;}
	
	/**
	 * Get the reduction in radial conductance due to ABA, environment and water potential
	 * @category water
	 * @category solute
	 * @return the reduction in radial conductance
	 */
	public float getLrPercent(){
		return (float) (getEnviModifier() * getStomataModifier()) ;
//		return (float) getABAModifier();
		}
	
	/**
	 * Return the photosynthesis over the last time step
	 * @category carbon
	 * @return the photosynthesis
	 */
	public double getPhotosynthesis(){return getPhotosynthesisFarquar();}
	
	
	/**
	 * Return the photosynthesis over the last time step based on the Farquar model
	 * V.max = maximum rubisco-limited rate in micromoles per (m^2 sec)
	 * J.max = maximum light-limited rate in micromoles per (m^2 sec)
	 * APAR = absorbed photosynthetically-active radiation in micromoles per (m^2 sec)
	 * c.i = intercellular CO2 partial pressure in Pascals (roughly ppm/10)
	 * @category carbon
	 * @return the photosynthesis
	 */
	public double getPhotosynthesisFarquar(){
		
		float APAR = (float) ExogenousEnvironment.getLight(node, parms.maxPAR);
		float Ci = 15 ;
		float stoMod = (float) Math.min(1, getStomataModifier()*2); // Stomatal effect on Carbon is less important
		
		Ci = Ci * stoMod;
		
		float timeSec = parms.timeStep * 3600;
		float surfaceM = (float) (this.getSurface() / 5e5);
			  
		// Some local parameters we need
		float psfc = 101325;  // surface air pressure (Pascals)
		float gamma = 3;  // CO2 compensation point (Pascals)
		float Oi = 0.209f * psfc;  // oxygen partial pressure in chloroplast
		float Kc = 30;  //Michaelis-Menten constant for carboxylation (Pascals)
		float Ko = 30000;  // Michaelis-Menten constant for oxidation (Pascals)
			  
		// Solution of quadratic (Bonan 17.8)
		float a = 0.7f;
		float b = -(Jmax + 0.385f * APAR);
		float c = 0.385f * Jmax * APAR;
		float J1 = (float) (-b + Math.sqrt(Math.pow(b, 2) - 4 * a * c) ) / (2 * a);
		float J2 = (float)(-b - Math.sqrt(Math.pow(b, 2) - 4 * a * c) ) / (2 * a);
		float J = Math.min(J1, J2);
			  
		// Rubisco-limited rate of photosynthesis
		float wc = Vmax * (Ci - gamma) / (Ci + Kc * (1 + Oi / Ko)); //  # Bonan 17.6
			  
		// Light-limited rate of photosynthesis
		float wj = J * (Ci - gamma) / (4 * (Ci + 2 * gamma)); //            # Bonan 17.7
			  
		// Sink-limited rate of photosynthesis
		float ws = Vmax / 2;
			  
		// Dark respiration 
		float Rd = 0.015f * Vmax;
			  
		// Net assimilation in µmol CO₂ m-2 s-1
		float An = Math.min(Math.min(wc, wj), ws) - Rd;
			  
		// Net carbon fixed
		float molMassCo2 = 44;  // molecular mass of CO₂
		float tot_carbon = (An / 1000) * timeSec * surfaceM  * molMassCo2;
			  
//		System.out.println("sto = "+this.getLrPercent());
		
		System.out.println(Math.max(tot_carbon * getEnviModifier(), 0));
		
		
		return Math.max(tot_carbon * getEnviModifier(), 0);
		
	}
	
	
	/**
	 * Update the photosynthesis value over the last time step
	 * It has to be done this way because water resolution can be done more than
	 * one time by time step.
	 * @category carbon
	 */
	public void updatePhotosynthesis(){
		double wue = Util.getInterpolatedValue(parms.wueX, parms.wueY, (float) getLrPercent());
		photosynthesis += wue * Math.abs(getRadialWaterFlux() * (parms.timeStep * 3600 /parms.resolveIteration) * 1e6);
	}
	
	/**
	 * Reset the photosynthesis value of this article to 0
	 * @category carbon
	 */
	public void resetPhotosynthesis(){ photosynthesis = 0;}

	/**
	 * Get the order of the article. For leaf, the leaf number
	 * @category architecture
	 * @return the leaf number
	 */
	public int getOrdre(){return this.leafNum;}

	/**
	 * Get the growth of the article
	 * @category architecture
	 * @return the growth in lenght
	 */
	public float getGrowth() {return (float) getAxialGrowth();}

	/**
	 * Get the distance from the base of this article
	 * For leaf, return 0
	 * @category architecture
	 * @return 0
	 */
	public float getDistFromBase() {return 0;}

	/**
	 * Get the aquaporin effect on the radial conductivity
	 * For leaf, no effect, return 1
	 * @category water
	 * @return 1
	 */
	public double getAQPCondEffect() {return 1;}
	
	/**
	 * Get the width of the segment
	 * @category architecture
	 * @return the width
	 */
	public double getWidth(){ return width;} 
  
	/**
	 * Get the length of the segment
	 * @category architecture
	 * @return the length
	 */
	public double getLength(){ return length;}

	/**
	 * Get the thickness of the segment
	 * @category architecture
	 * @return the thickness
	 */
    public double getThickness(){ return thickness;} 
    
	/**
	 * Get the surface of the segment
	 * @category architecture
	 * @return the surface
	 */
    public double getSurface(){ return width * length;}
  
	/**
	 * Get the volume of the segment
	 * @category architecture
	 * @return the volume
	 */
    public double getVolume(){ return width * length * thickness;} 	
	
}
