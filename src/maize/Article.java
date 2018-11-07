
package maize;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * This class defines the article, the building block of the plant network.
 * It is linked to both a netwrok and endogenous environment.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */
public abstract class Article {

	// General parameters
	public ParameterSet parms;			// PlanetCT parameter set	 
	public boolean isSelected = false;			// Flag for the selection of the article in CrossTalk
  	public double baseTemp; 					// Base temperature for the growth of the plant 

	// Topological parameters
	public Network network;        				// Network of the article (it is the base of the network)
  	public int internode = 0;					// Number of the internode (if relevant)
  	public int type = 0;						// Type of the article [10 = primary root, 20 = seminal root, 30 = crown root, 40 = stem, 41 = leaf]
  	public int ordre;							// Topological order of the article (1 = parent, 2 = first-order child, ...)
	
  	// Time parameters
	public double formationDate;     			// Formation date of the article [h]
	public double thermalAge;        			// Age of the article in thermal time [°Ch]
 	
	// Morphological parameters
	public EndogenousEnvironment envEndo;    	// Endogenous environment of the article
	public Point3d node;      					// Node located on the proximal part of the article (attached to the parent)
	public Vector3d axialdirection;     		// Axial direction of the article (dynamic, must be recalculate over time) 
	public Vector3d normalDirection;     		// Normal direction of the article (dynamic, must be recalculate over time) 
	public double relativePosition;       		// Relative position on parent article (0.0 proximal, 1.0 distal)
	public double insertionAngle;       		// Insertion angle, respectivelly to its carrier
	public double radialAngle;      			// Radial angle, respectivelly to its carrier
	public double spinAngle;   					// Spin angle of the article (rotation around its axis)
	public double diameter;						// Diameter of the article
  	public double length;						// Length of the article
	public double dryMass = 0;					// Dry mass of the articlr
  
	// Carbon parameters
	public double growthEfficiency = 1;					// Growth efficient coefficient depending on the carbon availability [0-1]
	public float carbonGrowthEfficiency = 0.73f;		// Amount of carbon needed for the growth of the article [gCO2/gDM] (Drouet & Pagès 2003)
	public float maintenanceRespirationRate = 0.0032f;	// Amount of carbon needed for the maintenance of the [gCO2/gDM]  (Drouet & Pagès 2003) 
	public int tempCoefficient = 2;						// Temperature coefficient used to compute the carbon maintenance demand  (Drouet & Pagès 2003)
	public int tRef = 20;								// Reference temperature used to compute the carbon maintenance demand [°C]  (Drouet & Pagès 2003)
	
	
	// Water parameters
  	public double KhEqu = 0.0 ;							// Equivalent resistance of the article
			 
  	// Space-time info parameters
  	static Soil soil = null;
  	static SpaceTimeInfo stiRH = null;
  	static SpaceTimeInfo stiWC = null;
  	static SpaceTimeInfo stiWP = null;
  	static SpaceTimeInfo stiHum = null;
  	static SpaceTimeInfo stiTemp = null;
  	static SpaceTimeInfo stiRM = null;


  	/**
  	 * Basic constructor
  	 */
  	public Article() {
  		
  		this.node=new Point3d();
  		this.axialdirection=new Vector3d();
  		this.normalDirection=new Vector3d();
  		this.relativePosition=0.0;
  		this.insertionAngle=0.0;
  		this.radialAngle=0.0;
  		this.spinAngle=0.0;
  		this.formationDate=Time.getTime();
  		this.thermalAge=0.0f;
  		this.envEndo=new EndogenousEnvironment();
  		this.network=null;
  		baseTemp = parms.baseTemp;
  	}
  
  	/**
  	 * Single constructor for the subsequent calculation of the positioning 
  	 * of the article used by the meristem to create metamere
  	 * @param posRel
  	 * @param angIns
  	 * @param angRad
  	 * @param angAzimuth
  	 * @param parms
  	 */
  	public Article(double posRel, double angIns, double angRad, double angAzimuth, ParameterSet parms) {
  		this.node=new Point3d();
  		this.axialdirection=new Vector3d();
  		this.normalDirection=new Vector3d();
  		this.relativePosition=posRel;
  		this.insertionAngle=angIns;
  		this.radialAngle=angRad;
  		this.spinAngle=angAzimuth;
  		this.formationDate=Time.getTime();
  		this.thermalAge=0.0f;
  		this.envEndo=new EndogenousEnvironment();
  		this.network=null;
  		this.parms=parms;
  		baseTemp = parms.baseTemp;
  	}
  
  	/**
  	 * 
  	 * @param posRel
  	 * @param angIns
  	 * @param angRad
  	 * @param angAzimuth
  	 * @param ordre of the article
  	 * @param type
  	 * @param parms
  	 */
  	public Article(double posRel, double angIns, double angRad, double angAzimuth, 
		  int ordre, int type, ParameterSet parms) {
  		this.node = new Point3d();
  		this.axialdirection = new Vector3d();
  		this.normalDirection = new Vector3d();
  		this.relativePosition = posRel;
  		this.insertionAngle  =  angIns;
  		this.radialAngle = angRad;
  		this.spinAngle = angAzimuth;
  		this.ordre = ordre;
  		this.type = type;
  		this.formationDate = Time.getTime();
  		this.thermalAge = 0.0f;
  		this.envEndo = new EndogenousEnvironment();
  		this.network = null;
  		this.parms = parms;
  		baseTemp = parms.baseTemp;
    }
  
  	/**
  	 * Constructor practice for propagules
  	 * Position, direction and endogenous environment known
  	 * @param node
  	 * @param dirAx
  	 * @param dirNorm
  	 * @param envEndo
  	 * @param parms
  	 */
  	public Article(Point3d node, Vector3d dirAx, Vector3d dirNorm, 
		  EndogenousEnvironment envEndo, ParameterSet parms) {
  		this.node = node;  
  		this.axialdirection = dirAx;
  		this.normalDirection = dirNorm;
  		this.relativePosition = 0.0;
  		this.insertionAngle = 0.0;
  		this.radialAngle = 0.0;
  		this.spinAngle = 0.0;
  		this.formationDate = Time.getTime();
  		this.thermalAge = 0.0f;
  		this.envEndo = envEndo;
  		this.network = null;
  		this.parms = parms;  
  		baseTemp = parms.baseTemp;
  	}


  	/**
  	 * Constructor
  	 * @param node
  	 * @param dirAx
  	 * @param dirNorm
  	 * @param envEndo
  	 * @param posRel
  	 * @param dateForm
  	 * @param parms
  	 */
  	public Article(Point3d node, Vector3d dirAx, Vector3d dirNorm, EndogenousEnvironment envEndo, 
		  double posRel, double dateForm, ParameterSet parms) {
  		this.node=node;
  		this.axialdirection=dirAx;
  		this.normalDirection=dirNorm;
  		this.envEndo=envEndo;
  		this.relativePosition=posRel;
  		this.formationDate=dateForm;
  		this.thermalAge=0.0f;
  		this.network=null;
  		this.parms=parms;
  		baseTemp = parms.baseTemp;
  	}
  
  
  	/**
  	 * Initialize the article, especially the space time info
  	 * @param s
  	 */
  	public static void initialize() {
  		Util.log("Initialization of planetMaize.Article");
  	}
  	
  	
//////////////////////////////////// GET FUNCTIONS /////////////////////////////////////////
  	
  	/**
  	 * Get the node associated with the article
   	 * @category architecture
 	 * @return the node
  	 */
  	public Point3d getNode() { return node; }

  	/**
  	 * Get the order of the article
  	 * @category architecture
  	 * @return the order
  	 */
  	public int getOrdre() { return ordre; }

  	/**
  	 * Get the axial direction of the article
  	 * @category architecture
  	 * @return the axial direction
  	 */
  	public Vector3d getAxialdirection() { return axialdirection;}

  	/**
  	 * Get the parent article of this one
  	 * @category architecture
  	 * @return the parent article
  	 */
  	public Article getParentArticle() {
  		if (network.parentNetwork!=null) return network.parentNetwork.baseArticle;
  		else return null;
  	}
  	
  	/**
  	 * Get the age, in hour, of the article
  	 * @category architecture
  	 * @return the age [h]
  	 */
  	public double getAge() { return (double) (Time.getTime()-formationDate); }
  	
  	/**
  	 * Get the increment in sum of thermal time based on the time step and the position of the article
  	 * Uses a trapezoidale rule to estimate the integral under the temp function
  	 * @category environment
  	 * @return the increment in thermal time
  	 */
    public double getDeltaST() {
    	double st;	
    	if(node.z > 0){
    		double cT = ExogenousEnvironment.getTemperature(Time.getCurrentHour() ,getCentralPosition(), stiTemp, parms.maxTemperature); 	// current temperature
    		double pT = ExogenousEnvironment.getTemperature(Time.getPreviousHour() ,getCentralPosition(), stiTemp, parms.maxTemperature);	//previous temperature
    		st = Time.getTimeStep()*((cT+pT)/2);
    		return st;		
    	} 
    	else return Time.getTimeStep() * (Soil.getTemperature(node)-baseTemp);
    }  
   
  	/**
  	 * Get the central position of the article as a Point3D
  	 * @category architecture
  	 * @return the center of the article
  	 */
  	public Point3d getCentralPosition() {
  		Point3d centre=new Point3d(axialdirection);
  		centre.scaleAdd(getLength() / 2.0, node);
  		return centre;
  	}
    
  	/**
  	 * Get the cavitation effect on the axial conductance of the article based on is water potential
	 * Weibull function between 0 and 1
	 * From Li et al 2009 / Sperry et al 2002
  	 * @category water
  	 * @return the cavitation effect
  	 */
	public double getCavitationEffect(){ 
		
		float[] parm = null;
		boolean flag = true;
		
		// Set the parameters depending on the article type
		if(type == 10 || type == 20 || type == 30) parm = new float[]{2.8f, 1.4f};		
		else if(type == 11 || type == 21 || type == 31) parm = new float[]{2.8f, 0.9f};		
		else if(type == 40 || type == 41) parm = new float[]{3f, 1.8f};
		else flag = false;
		
		double x = -getWaterPot() - parms.cavitationMod;
		if(x < 0) x = 0;
		
		if(parms.cavitation && flag) {
//			Util.log("cavitation = "+Math.exp(-Math.pow(((-getWaterPot())/parm[1]),parm[0])));
			float cavitation = (float) Math.exp(-Math.pow((x/parm[1]),parm[0]));
			if(cavitation == 0) return Util.epsilonFloat;
			else return cavitation;
		}
		else return 1;
		
	}
	
	/**
	 * Get the diameter of the article
  	 * @category architecture
	 * @return the diameter
	 */
	public double getDiameter(){return diameter;}
	
	/**
	 * Get the length of the article
  	 * @category architecture
	 * @return the length
	 */
	public double getLenght(){return length;}
  	
	/**
	 * Get the growth efficient coefficient of the article.
	 * This coefficient is defined by the C availability
	 * @category carbon
	 * @return the growth efficiency coefficient
	 */
	public double getGrowthEfficiency(){return growthEfficiency;}
	 
	/**
	 * Get the difference in a given solute quantity [nMol] between
	 * this step and the previous one
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flag [1 = ABA]
	 * @return the difference in the solute quantity 
	 */
	public double getDeltaQuantSolute(int flag) {
		return envEndo.getDeltaQuantSolute(flag);
		}

	/**
	 * Get how much solute has been displaced to the father element.
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flag [1 = ABA]
	 * @return the flow of solute
	 */
	public double getFlowSolute(int flag) {
		return envEndo.getQuantFlowSolute(flag);
	}
	
	/**
	 * Get the difference in a given solute concentration between
	 * this step and the previous one
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flag [1 = ABA]
	 * @return the difference in the solute concentration 
	 */
	public double getDeltaConcSolute(int flag){
		return envEndo.getDeltaConcSolute(flag);
	}
	
	/**
	 * Get the concentration of a given solute (in nM), function of quantSolute and the volume.
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flag [1 = ABA]
	 * @return the concentration
	 */
	public double getConcSolute(int flag) {
		return envEndo.getQuantSolute(flag) / (getVolume() * 10.0e-3);
	}
	
	/**
	 * Get the mass of a given solute [g].
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flag [1 = ABA]
	 * @return the mass
	 */
	public double getMassSolute(int flag) {
		return getVolume()*getConcSolute(flag)*ParameterSet.molMassABA*10e-3;
	}	
	
	/**
	 * Get the water potential [MPa] of the article
	 * @category water
	 * @return the water potentiel
	 */
	public double getWaterPot() { return envEndo.getWaterPot(); }
	
	/**
	 * Get the water radial resistivity (called here resistance) of the article [MPa s m-1].
	 * Resistance = 1 / conductance 
	 * @category water
	 * @return the water radial resistance
	 */
	public double getWaterRadialResistance(){ return 1 / getWaterRadialConductance(); };
	
	/**
	 * Get the water radial resistivity (called here resistance) of the article [MPa s m-1].
	 * Resistance = 1 / conductance 
	 * @category water
	 * @return the water radial resistance
	 */	
	public double getWaterAxialResistance(){ return 1 / getWaterAxialConductance(); }
	
	/**
	 * Get the radial water flux in the article
	 * @category water
	 * @return the water radial flux
	 */
	public double getRadialWaterFlux(){ return envEndo.getRadialWaterFlux(); }
	  
	/**
	 * Get the radial water flux in the article
	 * @category water
	 * @return the water radial flux
	 */ 
	public double getAxialWaterFlux(){ return envEndo.getAxialWaterFlux(); }
	    
	/**
	 * Get the water potential outside of the article.
	 * NOTE: sti does not seems to be able to be passed...
	 * @category water
	 * @return the exogenous water potentiel
	 */
	public double getWaterPotExo(){
//		if(parms.stressType == 0 && node.z < 0) return soil.getDouble(stiWP.setCoordinates(node.x, node.y, node.z));
//		else 
		return ExogenousEnvironment.getWaterPotExo(node, parms.stressType, parms.startStress, parms.reWateringFrequence);
	}
	  
	/**
	 * Get the type of the article, which define who it is
	 * [10 = primary root, 20 = seminal root, 30 = crown root, 40 = stem, 41 = leaf]
	 * @category architecture
	 * @return type
	 */
	public int getType() {
		return type;
	}
	
	/**
	 * Get the dry mass of the article
	 * @category carbon
	 * @return the dry mass of the article
	 */
	public double getDryMass(){return dryMass;}
								
	/**
	 * Check if this article is selected in CrossTalk
	 * @return true if selected
	 */
	public boolean isSelected(){
		return isSelected;
	}
	
	////////////////////////////////////  SET FUNCTIONS  /////////////////////////////////////////
		
	/**
	 * Set if this article is selected in CrossTalk
	 * @param s true is selected
	 */
	public void setSelected(boolean s){
		isSelected = s;
	}
		
	/**
	 * Set the quantity of a given solute based on its concentration
	 * @author Vincent Larondelle
	 * @category solute
	 * @param conc the concentration
	 * @param flag [1 = ABA]
	 */
	public void setQuantSolute(double conc,int flag) {
		double quant = conc * getVolume() * 10.0e-3;
		envEndo.setQuantSolute(quant,flag);
	}
	
	/**
	 * Set the quantitative flow of a solute in direction of the father (negative value
	 * means from the father in direction of this article)
	 * @author Vincent Larondelle
	 * @category solute
	 * @param flow
	 * @param flag
	 */
	public void setQuantFlowSolute(double flow,int flag) {
		if (flag == 1) this.envEndo.setQuantFlowSolute (flow,flag);
	}

	/**
	 * Set the difference in quantity of a given solute based on its new anf former values
	 * @author Vincent Larondelle
	 * @category solute
	 * @param quant
	 * @param newQuant
	 * @param flag [1 = ABA]
	 */
	public void setDeltaQuantSolute(double quant, double newQuant, int flag) {
		envEndo.setDeltaQuantSolute(newQuant - quant, flag);
	}
	
	/**
	 * Set the difference in concentration of a given solute based on its new anf former values
	 * @author Vincent Larondelle
	 * @category solute
	 * @param conc
	 * @param newConc
	 * @param flag [1 = ABA]
	 */
	public void setDeltaConcSolute(double conc, double newConc, int flag){
		envEndo.setDeltaConcSolute(newConc - conc, flag);
	}
	
	/**
	 * Set the water potential of the article
	 * @category water
	 * @param pot the water potential
	 */
	public void setWaterPot(double pot){
		envEndo.setWaterPot(pot);
	}
	
	/**
	 * Set the axial water flux in the plant as defined by the Solver [m3/s]
	 * @category water
	 */
	public void setAxialWaterFlux(){ 
		double flux = 0;
			
		if(getParentArticle() != null){
			flux = -1 * (getParentArticle().envEndo.getWaterPot() - envEndo.getWaterPot()) * (1 / getWaterAxialResistance()) / (getLength() / 1.0e2);
		}
		// the virtual article has no father, hence the flux is defined by his children
		else{
			if (network.childNetworkList != null) {
				for (Network resFils : network.childNetworkList) {
					if(resFils.baseArticle instanceof SegStem){
						flux = -1 * (envEndo.getWaterPot()-resFils.baseArticle.envEndo.getWaterPot()) * (1 / resFils.baseArticle.getWaterAxialResistance()) / (resFils.baseArticle.getLength() / 1.0e2);
					}
				}
			}
		}
		envEndo.setAxialWaterFlux(flux);
	};
	  
	/**
	 * Set the axial water flux of the article
	 * @category water
	 * @param axialFlux
	 */
	public void setAxialWaterFlux(double axialFlux){
		envEndo.setAxialWaterFlux(axialFlux);
	};   
	 
	/**
	 * Set the radial water flux in the plant as defined by the Solver [m3/s]
	 * @category water
	 * @param axialFlux
	 */
	public void setRadialWaterFlux(){
		envEndo.setRadialWaterFlux((1.0/getWaterRadialResistance())*(getWaterPotExo()-getWaterPot())*(getSurface()/1.0e4));
	};
	
	/**
	 * Set the radial water flux of the article
	 * @category water
	 * @param axialFlux
	 */
	public void setRadialWaterFlux(double radialFlux){
		envEndo.setRadialWaterFlux(radialFlux);
	};
	  
	/**
	 * Set the growth efficiency coefficient of the article
	 * @category carbon
	 * @param g: the growth efficiency coefficient
	 */
	public void setGrowthEff(double g){growthEfficiency = g;}
	
  	/**
  	 * Set the node of the article
  	 * @category architecture
  	 * @param nd
  	 */
  	public void setNode(Point3d nd) { node=nd;}

  	/**
  	 * Set the axial direction of the article
  	 * @category architecture
  	 * @param dirAx
  	 */
  	public void setAxialDirection(Vector3d dirAx) { this.axialdirection=dirAx; }

  	/**
  	 * Set the formation of the article
  	 * @category architecture
  	 * @param dateForm
  	 */
  	public void setDateForm(double dateForm) { this.formationDate=dateForm; }

  	/**
  	 * Set the relative position of the article
  	 * @category architecture
  	 * @param posRel
  	 */
  	public void setRelativePosition(double posRel) { this.relativePosition=posRel;}

  	/**
  	 * Set the network of the article
  	 * @category architecture
  	 * @param net
  	 */
  	public void setNetwork(Network net) { this.network=net; }
  	
  	/**
  	 * Rotate the article around its axis
   	 * @category architecture
 	 * @param angRot
  	 */
  	public void applySpinRotation(double angRot) {
  		spinAngle+=angRot;
  		while (spinAngle>(2.0*Math.PI)) spinAngle-=2.0*Math.PI;
  	}

  	/**
  	 * Apply a deviation in x, y, z on the article
  	 * @category architecture
  	 * @param x
  	 * @param y
  	 * @param z
  	 */
  	public void applyDeviation(double x, double y, double z) {applyDeviation(new Vector3d(x, y, z));}
  
  	/**
  	 * Apply a deviation in the given direction on the article
   	 * @category architecture
 	 * @param direction
  	 */
  	public void applyDeviation(Vector3d direction) {
  		Vector3d newAxialDir = new Vector3d(axialdirection);
  		newAxialDir.add(direction);
  		newAxialDir.normalize();

  		Article parent = network.parentNetwork.baseArticle;
  		insertionAngle = parent.axialdirection.angle(newAxialDir);
     
  		Vector3d tmp = new Vector3d(); 
  		tmp.scaleAdd(-newAxialDir.dot(parent.axialdirection), parent.axialdirection, newAxialDir);

  		Vector3d tmp2 = new Vector3d();
  		tmp2.cross(parent.normalDirection, tmp);
  		double sign = Math.signum(tmp2.dot(parent.axialdirection));
  		radialAngle = sign * parent.normalDirection.angle(tmp);
  	}
  

  	/**
  	 * Compute the orientation (axial and normal direction) of the article
  	 * based on its parent and using the angles radialAngle, insertionAngle and spinAngle
  	 * @category architecture
  	 * @param parent
  	 */
  	public void setOrientation(Article parent) {
    
  		if (parent!=null) {

  			Vector3d parentAxDir = parent.axialdirection; 		// Parent axial direction, A	 	
  			Vector3d parentNormDir = parent.normalDirection;  	// Parent normal direction, N
  			
  			// Rotation of N around A of the radial angle -> nR
  			Vector3d nR = new Vector3d(parentNormDir); 
  			Util.rotation3D(nR, parentAxDir, radialAngle);
  			
  			// Calculation of the rotation axis
  			Vector3d oO = new Vector3d();
  			oO.cross(parentAxDir, nR);
  			
  			// Rotation of A around oO of the insertion angle(-> A')
  			axialdirection.set(parentAxDir);
  			Util.rotation3D(axialdirection, oO, insertionAngle);

  			if (ordre == parent.ordre) {
  				// Rotation of N around oO of the insertion angle (-> N')
  				normalDirection.set(parentNormDir);
  				Util.rotation3D(normalDirection, oO, insertionAngle);
  			} 
  			else {
  				// Positionning of N' in a plane formed by A and A', oriented through A
  				// If A and A' are colinear, keep the normal direction of the parent (normally should not append)
  				double mu = axialdirection.dot(parent.axialdirection);
  				if (Math.abs(mu) > 0.9999) normalDirection.set(parent.normalDirection);
  				else {
  					normalDirection.scaleAdd(-mu, axialdirection, parent.axialdirection);
  					normalDirection.normalize();
  				}
  			}
      
  			// Test if N' is in the plane defined by A and A' if it is a ramification
  			if (ordre != parent.ordre) {
  				Vector3d checkV = new Vector3d();
  				checkV.cross(parent.axialdirection, axialdirection);
  				if (Math.abs(normalDirection.dot(checkV)) > 1e-10)
  					Util.log("Ramification alignement error: ordres " + ordre + " on " + parent.ordre+"");
  			}
 
  			// Rotation of N' around A' of spin angle
  			Util.rotation3D(normalDirection, axialdirection, spinAngle);
  		}
  		else Util.log("Problem: parent == null (Article.setOrientation)"); 		
  	}

  	/**
  	 * Calculates and assigns the position of the article according
  	 * to the parent position and direction
  	 * @category architecture
  	 * @param parent
  	 */
  	public void setPosition(Article parent) {
  		if (parent!=null) {
  			setOrientation(parent);
  			node.scaleAdd(relativePosition * parent.getLength(), parent.axialdirection, parent.node);
  		}
  		else Util.log("Problem: parent == null (Article.positionne)");
  	}

 
	////////////////////////////////////  UPDATE FUNCTIONS  /////////////////////////////////////////
	
    /**
     * 	Update the thermal age of the article
   	 * @category architecture
     */
    public void updateThermalAge() { thermalAge += getDeltaST(); }
  	
  
    /**
     * Computation of the reaction component of the variation of solute quantity in an article
     * @author Vincent Larondelle
     * @category solute
     * @param flag [1 = ABA]
     * @return
     */
    public double soluteReactBalance(int flag){
    	if (flag==1) return getABAProduction() - getABAConsumption();
    	else return 0;
    }
  
    /**
     * Update quantity of a given solute based on its concentration and a coeff factor
     * @author Vincent Larondelle
     * @category solute
     * @param coeff
     * @param flag [1 = ABA]
     */
    public void updateSoluteQuant(double coeff, int flag) {
    	setQuantSolute(getConcSolute(flag) * (1.0 + coeff), flag);
    }
    
    /**
     * Modify the water content of the exogenous environment
     * @category water
     * @param WC
     */
    public void updateWaterContExo(double WC){
    	// If hydrolic lift is not allowed, the plant can not rewater the soil
//    	if(!parms.allowHydraulicLift && WC < 0) return;
//    	if(node.z < 0) {
//    		if(parms.splitRootType == 0) soil.addDouble(stiWC.setCoordinates(node.x, node.y, node.z), -WC);
//    		if(parms.splitRootType == 1) if(node.x > parms.splitRootLimit) soil.addDouble(stiWC.setCoordinates(node.x, node.y, node.z), -WC);
//    		if(parms.splitRootType == 2) if(node.z > -parms.splitRootLimit) soil.addDouble(stiWC.setCoordinates(node.x, node.y, node.z), -WC);
//    	}  
    }
	
	/**
	 * Initialize the concentration of the "new" quantity of a solute
	 * @author Vincent Larondelle
	 * @category solute 
	 * @param flag [1 = ABA]
	 */
	public void initSoluteNew(int flag) {
		envEndo.initSoluteNew(flag);
	}
	
	/**
	 * Initialises the "quantSoluteNew", gets the reaction balance and update the "new" solute value
	 * @author Vincent Larondelle
	 * @category solute [1 = ABA]
	 * @param flag [1 = ABA]
	 */
	public void soluteReaction(int flag) {
		initSoluteNew(flag);
		envEndo.setQuantSolute(envEndo.getQuantSolute(flag+1) + soluteReactBalance(flag), flag);
//		Util.log("react = "+soluteReactBalance(flag));
//		Util.log("new = "+(envEndo.getQuantSolute(flag+1) +soluteReactBalance(flag)));
//		Util.log("-----------------");
	}  
	
	/**
	 * Kill this article and all its children. They will be removed from the network
	 */
	public void killArticle(){
		int num = this.network.getNumberOfChildren();
		if(this.network.childNetworkList != null){
			this.network.removeAllChildNetwork();
			Util.log("Excission successfull: "+num+" article(s) removed");
		}
		else Util.log("No article to remove");
	}

	////////////////////////////////////  ABSTRACT FUNCTIONS  /////////////////////////////////////////
	
	/**
	 * Get the volume, depend on the geometry of the article
	 * @category architecture
	 * @return the volume [cm3]
	 */
	public abstract double getVolume();
	
	/**
	 * Get the length of the article
	 * @category architecture
	 * @return the length [cm]
	 */
	public abstract double getLength(); 
	
	/**
	 * Get the width of the article
	 * @category architecture
	 * @return the width [cm]
	 */
	public abstract double getWidth();
	
	/**
	 * Get the thickness of the article
	 * @category architecture
	 * @return the thickness [cm]
	 */
	public abstract double getThickness();
	
	/**
	 * Get the surface of the article
	 * @category architecture
	 * @return the surface [cm2]
	 */
	public abstract double getSurface();
	
	/**
	 * Get the length of the article
	 * @category architecture
	 * @return the length [cm/°Ch]
	 */
	public abstract float getGrowth(); 

	/**
	 * Get the distance from the base of the network 
	 * @category architecture
	 * @return the distance to the base [cm]
	 */
	public abstract float getDistFromBase();
	
	/**
	 * Get the article identifier. Depend on the article type
	 * @category architecture
	 * @return the article identifier
	 */
	public abstract int getArticleID(); 

	/**
	 * Insure the devellopement of the plant (age, size, mass, ...)
	 * @category architecture
	 */
	public abstract void develop(); 
	    	    
	/**
	 * Get the ABA produced by the article during the time step
	 * @category solute
	 * @return the ABA production
	 */
	public abstract double getABAProduction();	  
	
	/**
	 * Get the ABA consumption of the article
	 * @category solute
	 * @return the ABA condumption
	 */
	public abstract double getABAConsumption();
	  
	/**
	 * Get the effect of AQP on the radial conductance.
	 * Only for roots, else returns 0
	 * @category water
	 * @return the AQP effect on the radial conductance
	 */
	public abstract double getAQPCondEffect();
	  		
	/**
	 * Get the reduction in root radial conductance.
	 * Only for leaf
	 * @category water
	 * @return the reduction in radial conductance
	 */
	public abstract float getLrPercent();
				
	/**
	 * Get the equivalent conductance of the article
	 * @category architecture
	 * @return the equivalent resistance
	 */
	public abstract double getEquivalentConductance();
	
	/**
	 * Get the water radial conductance of the article
	 * @category water
	 * @return the water radial conductance
	 */
	public abstract double getWaterRadialConductance(); //return 1.0e-10;
			
	/**
	 * Get the water axial conductance of the article
	 * @category water
	 * @return the water axial conductance
	 */
	public abstract double getWaterAxialConductance();//{return 2.5e-8;//2.5e-8  

	/**
	 * Get the maintenance demand in C of the article
	 * @category carbon
	 * @return the amount of C required for maintenance [gC]
	 */
	public abstract double getMaintenanceDemand();

	/**
	 * Get the growth demand in C of the article
	 * @category carbon
	 * @return the amount of C required for growth [gC]
	 */
	public abstract double getGrowthDemand();
	
	/**
	 * Get the photosynthesis production in c for the article
	 * @category carbon
	 * @return the amount of C produce by photosynthesis [gC]
	 */
	public abstract double getPhotosynthesis();

  }

///**
// * Get the maintenace demand of the article
// * From Drouet and Pagès 2003
// * @return
// */
//public abstract double getMaintenanceDemand(){
//	float temp = (float) (Environment.getTemperature(node, stiTemp) - tRef) / 10;	// Drouet and Pages
//	return maintenanceRespirationRate * dryMass * Math.pow(tempCoefficient, temp) * getDeltaST();	// Drouet and Pages
//}


