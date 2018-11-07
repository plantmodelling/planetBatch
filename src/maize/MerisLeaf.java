package maize;


/**
 * This class represent a leaf meristem.
 * Its aims is to generate new leaf segments.
 * Is located between the sheath and the lamina, at the ligula
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium) 
 */

public class MerisLeaf extends Article {
  
	public double growthEndAgeLamina;  		// Limit age (thermal age) for the growth of the lamina
	public double growthEndAgeSheath;		// Limit age (thermal age) for the growth of the sheath

	// Carbon
	public double leafPotGrowth ;			// Potential growth of the leaf meristem
	public double dryMassRatio = 0.001f;	// Dry mass ratio (gDM/volume)
  	
	// ABA
	public float ABAHalfLife = 0.96f;		// Half life of ABA [h] (Ren et al. 2007)
	
	// Architecture
	public double volInit;	    			// Initial volume of the article
	public double diameter;       			// Diameter of the meristem 
	public double width;					// Width of the leaf
	public double widthSheath;				// Width of the Sheath
	public double totalLength;				// Total length of the leaf
	public double length = 0.01;			// Length of the meristem
	public double thickness = 0.01;				// Thickness of the leaf
	public float curve = 30;				// Curvature of the leaf. Just for display, not used forlight interseption
	public int leafNum;						// Identifier of the leaf  
	public boolean lamina = false;			// Is the lamina growing ? 
	public SegStem stem;					// The stem segment to which this leaf is attached
		
	// Leaf Growth. From Drouet and Pagès 2003
	public static float leafGrowthRate ; 

	public double fLengthLamina;			// final length of the lamina
	public float ala = -6.36f;				// shape parameter
	public float bla = -1.4f;				// shape parameter
	public float nM = 11.54f;				// leaf number of the longest leaf
	public int nT = 17;						// total number of leaf
	public int LiM = 78;					// length of the longest leaf
  
	public double fLengthSheath;			// final lenght of the sheath
	public double als = 0.8f;
	public double bls = 1.3f;
	public double cls = -2.1f;
	public double wiM = 11;
  
	public double aLs = 3.077;
	public double bLs = 20.48;
	public double cLs = 0.569;
	public int nLs = 6;
  
	public double DIE;
	public double adie = -5.16;
	public double bdie = 58.2;
	public double ndie = 3.65;
	public double cdie = 462.97;
	public double ddie = 15.64;
  

  
	/**
	 * Constructor
	 */
	public MerisLeaf() {
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
	 * @param l
	 * @param s
	 */
	public MerisLeaf(double posRel, double angIns, double angRad, double angSpin,
			float diam, int ordre, int type, ParameterSet parms, int l, SegStem s) {
		super(posRel, angIns, angRad, angSpin, ordre, type, parms);
		stem = s;
		leafNum = l;
		setParams();
		volInit = widthSheath * length * thickness;
	}
  
    
	/**
	 * Set the parameters for the meristem
	 * Based on informations from Girardin 98 and Drouet 03
	 * @category architecture
	 */
	public void setParams(){
	  
		leafGrowthRate = parms.leafGrowthRate;
	  
		if(leafNum < nM) DIE = adie + (bdie * leafNum);
		else DIE =cdie + (ddie * leafNum);
	  
		fLengthLamina = (float) (LiM * Math.exp(((ala / 2) * 
				Math.pow((leafNum/nM) - 1, 2)) + ((bla / 2) *
						Math.pow((leafNum/nM) -1, 3))));	
	  
		growthEndAgeLamina = fLengthLamina / leafGrowthRate;
	  
		setWidth();
	  
		if(leafNum <= nLs) fLengthSheath = (float) aLs * leafNum;
		else fLengthSheath = bLs - (cLs * leafNum);
		
		growthEndAgeSheath = DIE+(fLengthSheath / (float) leafGrowthRate );
		widthSheath = (float) width;	  
	}
  
	/**
	 * Set the width of the leaf segment based on its position from the base of the leaf and the leaf number
	 * From Drouet 2003
	 * @category architecture
	 */  
	public void setWidth(){
	  	wiM = (-0.02 * Math.pow(leafNum, 3)) + (0.4304 * Math.pow(leafNum, 2)) - (1.666 * leafNum) + 3.08 ; 
	  	width = wiM * (als + (bls * (1-(totalLength / fLengthLamina))) + (cls * Math.pow((1-(totalLength / fLengthLamina)), 2)));
  	}
  
	/**
	 * Insure the developement of the meristem by modifiyng its attributes based on its age, ...
	 * Create new article (SegLeaf) and anchor them on the global network.
	 * Active only during daylight
	 */
	public void develop() {
		
		if(Util.isDayTime()){		// If day time
			
			if(thermalAge > DIE && thermalAge < (growthEndAgeLamina + growthEndAgeSheath) ){  

  				Network n =  createLeafSegment();						// create the leaf segment
  				if(n != null) network.addProximalNetwork(n);	// Add the newly create segment to the network			
  					
  				// Update parameters
   				relativePosition = 1; 	      
  				spinAngle = 0.0;
  				insertionAngle = 0.0;
  				radialAngle = 0.0;
  				setWidth();
  			}	
   		    updateThermalAge();
		}
	}  	 
  	  	
	/**
	 * Create a newly leaf segment, which will be, depending on the age of the meristem, 
	 * part of the sheath or the lamina
	 * @category architecture
	 * @return the new leaf segment
	 */
	public Network createLeafSegment() {
  		  
		Article en;			
		// Sheath
		if(thermalAge <= growthEndAgeSheath && growthEfficiency > 0){
			en = new SegLeaf(relativePosition, insertionAngle, radialAngle, spinAngle, widthSheath, 0.01 * growthEfficiency, 
					0.01, ordre, type, parms, leafNum, lamina, leafGrowthRate * getDeltaST(), stem);
			return new Network(en);
		}
  			
		// Lamina
		else if(thermalAge > growthEndAgeSheath && growthEfficiency > 0){
  				
  			// The insertion angle of the newly created segment (for the rendering of a curved leaf)
			insertionAngle = (leafGrowthRate * getDeltaST()) / curve;
			lamina = true;

			en = new SegLeaf(relativePosition, insertionAngle, radialAngle, spinAngle, width, 0.01 * growthEfficiency, 
					0.01, ordre, type, parms, leafNum, lamina, leafGrowthRate * getDeltaST(), stem);
			totalLength += leafGrowthRate * getDeltaST();   
			return new Network(en);
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
	 * Get the quantitative production of ABA, based on the volume and the time step
	 * From Dodd 2008
	 * @author Vincent Larondelle
	 * @category solute
	 * @return the aba production
	 */
	public double getABAProduction() {
		return 1 * (-170 * getWaterPot()) * Time.getTimeStep() * getVolume(); 
	}
 	
  	
	/**
	 * Get the growth demand for the article at this time step
	 * Used to compute the over plant requirement for growth
	 * From Drouet and Pagès 2003
	 * @category carbon
	 * @return the growth demand
	 */
	public double getGrowthDemand(){
		
		if(Util.isDayTime()){
			if(thermalAge <= growthEndAgeSheath) return 0.01 * widthSheath * 0.01 * carbonGrowthEfficiency * dryMassRatio;
			else return 0.01 * width * 0.01 * carbonGrowthEfficiency * dryMassRatio;
		}
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
	 * Get the article diameter, in this case, the width of the sheath
	 * @category architecture
	 * @return the width of the sheath
	 */
	public double getDiameter() { return widthSheath; } 
  
	/**
	 * Get the article length
	 * @category architecture
	 * @return the length
	 */
	public double getLength() { return length; } 
  
	/**
	 * Get the article width, in this case, the width of the sheath
	 * @category architecture
	 * @return the width of the sheath
	 */
	public double getWidth() { return widthSheath; } 

	/**
	 * Get the article thickness
	 * @category architecture
	 * @return the thickness of the article
	 */
	public double getThickness() { return thickness;}
  
	/**
	 * Get the article volume, in this case, the volume of a 3D rectangle
	 * @category architecture
	 * @return the volume
	 */
	public double getVolume() { return widthSheath * length * thickness; } 

	/**
	 * Get the radial conductance to water flux of the segment.
	 * The conductance is low, the meristem is assumed not to
	 * take part in the transpiration.
	 * @category water
	 * @return the water radial conductance
	 */
	public double getWaterRadialConductance() {
		return Util.epsilonDouble;
	}

	/**
	 * Get the axial conductance of the article.
	 * In this case, related to the conductance of the associated stem segment
	 * @category water
	 * @return the water axial conductance
	 */
	public double getWaterAxialConductance() {
    	return (stem.getWaterAxialConductance() / 10 ) * getCavitationEffect(); 
//    	return (1 / getWaterAxialResistance()) * getCavitationEffect(); 
	}
	
	/**
	 * Get the axial resistance of the article.
	 * In this case, related to distance from the leaf base [Wei et al. 1999]
	 * @category water
	 * @return the water axial resistance
	 */
	public double getWaterAxialResistance() {
    	return (stem.getWaterAxialConductance() / 10 ) * getCavitationEffect(); 
	}
	
	/**
	 * Get the surface of the article
	 * @category architecture
	 * @return the surface
	 */
	public double getSurface(){ return widthSheath * length * 4; }

	/**
	 * Ge the equivalent resistance of the article, in the case, 0
	 * @category water
	 * @return 0
	 */
	public double getEquivalentConductance() { return 0; }

	/**
	 * Get the growth of the article, in this case, 0
	 * @category architecture
	 * @return 0
	 */
	public float getGrowth() { return 0;}

	/**
	 * Get the distance from the base of the network. 
	 * In this case, 0
	 * @category architecture
	 * @return 0	 
	 */
	public float getDistFromBase(){ return 0; }

	/**
	 * Get the article ID 
	 * In this case, 0
	 * @category architecture
	 * @return 0	 
	 */
	public int getArticleID() { return 0; }

	/**
	 * Get the article aquaporin effect 
	 * In this case, 1
	 * @category water
	 * @return 1	 
	 */
	public double getAQPCondEffect() { return 1; }

	/**
	 * Get the article Lr reduction 
	 * In this case, 1
	 * @category water
	 * @return 1	 
	 */
	public float getLrPercent() { return 1; }

	/**
	 * Get the article photosynthesis 
	 * In this case, 0
	 * @category carbon
	 * @return 0	 
	 */
	public double getPhotosynthesis() { return 0;}

}
