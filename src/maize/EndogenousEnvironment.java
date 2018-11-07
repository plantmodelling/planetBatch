package maize;


/**
 * Represents the endogenous environment associated with each article, 
 * which will help determine its development.
 * Encapsulates the set of variables that characterize this environment, 
 * ie the concentration of various substances (resources, signals) and the 
 * mechanical stresses (weight)
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public class EndogenousEnvironment {
  
  
	/* ABA methods
	 * 
	 * Many methods of this model work with ABA, but could also be used for other hormones and solutes. Therefore, these methods
	 * use the argument "flag", which is the reference of the solute the method has to work on.
	 * 	Flag 1 is used for ABA
	 * 	Flag 2 is used for ABANew, which means ABA on the next time step
	 * I recommend to use impair numbers for solutes, and the next pair numbers for their "new" values
	 * For example, if flag 3 is used for auxin, flag+1 = flag 4 would be used for auxinNew */
  
	private static double maxAdmisQuantABA = 1.0e3; 	// It is not allowed to have more than 1 µmol ABA in an article
	private static double minAdmisQuantABA = 0.0; 		// A negative value for quantABA is not allowed
	private double quantABA = minAdmisQuantABA; 		// This is the main variable for ABA. It is measured in nmol.
	private double quantABANew; 						// quantABA at the next time step.  
	private double flowABA; 							// How much ABA is displaced out of this article. nmol/heure 
	private double deltaConcABA, deltaQuantABA; 
	static final double minAdmisConcABA = 1.0;
	static final double maxAdmisConcABA = 10.0e10;		//If we consider a very small amount for an item, for example, the volume of a meristem, 0.1mm ³, 
														//and the greatest amount of ABA, 1000 nmol, was a concentration of 10 ^ 10 nM, ie 10 M. 

	// Water methods
	private static double waterPotMinAdmis=-0.4;
	private double waterPot = waterPotMinAdmis; 		// Initial water potential value of the article
	private double radialWaterFlux = 0.0; 				// [m³/s]
	private double axialWaterFlux = 0.0; 				// [m³/s]
  
	
	/**
	 * Constructor
	 */
	public EndogenousEnvironment() {}

	/**
	 * Initialize solute quantity at the specified value
	 * @param quant: the value
	 * @param flag:	the solute type
	 */
	public EndogenousEnvironment(double quant,int flag) {
		setQuantSolute(quant,flag);
	}


	/**
	 * Set the quantitative flow for the given solute
	 * @param flow
	 * @param flag [1 = ABA]
	 */
	public void setQuantFlowSolute(double flow,int flag) { 
		if (flag == 1) this.flowABA=flow; 
	}

	/**
	 * Set the difference in concentration for a given solute
	 * @param conc
	 * @param flag[1 = ABA]
	 */
	public void setDeltaConcSolute(double conc, int flag){
		if(flag == 1) this.deltaConcABA = conc;
	}
	
	/**
	 * Set the difference in quantity [nMol] for a given solute
	 * @param quant
	 * @param flag [1 = ABA]
	 */
	public void setDeltaQuantSolute(double quant, int flag){
		if(flag == 1) this.deltaQuantABA = quant;
	}

	/**
	 * Getthe difference in concentration for a given solute
	 * @param flag [1 = ABA]
	 * @return the difference in concentration
	 */
	public double getDeltaConcSolute(int flag){
		if(flag == 1) return deltaConcABA;
		else return 0;
	}
	
	/**
	 * Get the difference in quantity [nMol] for a given solute
	 * @param flag [1 = ABA]
	 * @return the difference in quantity
	 */
	public double getDeltaQuantSolute(int flag){
		if(flag == 1) return deltaQuantABA;
		else return 0;
	}
	
	/**
	 * Get the quantitative flow for a given solute
	 * @param flag [1 = ABA]
	 * @return the flow
	 */
	public double getQuantFlowSolute(int flag) {
		if (flag ==1) return flowABA; 
		else return 0;
	}

	/**
	 * Get the quantity of a given solute
	 * @param flag [1 = ABA, 2 = ABA new]
	 * @return the quantity
	 */
	public double getQuantSolute(int flag) {
		if (flag==1) return quantABA;
		if (flag==2) return quantABANew;
		else return 0;
	}

	/**
	 * Set the quantity of a given solute in the article
	 * @param quant
	 * @param flag  [1 = ABA, 2 = ABA new]
	 */
	public void setQuantSolute(double quant, int flag) {
		// ABA
		if (flag==1){ 
			if (quant < minAdmisQuantABA) {
				this.quantABA = minAdmisQuantABA;
				Util.log("ABA lower bound exceeded: deficit = "+(minAdmisQuantABA-quant));
			}
			else if(quant > maxAdmisQuantABA) {
				this.quantABA = maxAdmisQuantABA;
				Util.log("ABA upper bound exceeded: excess = "+(quant-maxAdmisQuantABA));
			}
			else this.quantABA = quant;
		}
		// ABA new
		if (flag==2){ 
			if (quant < minAdmisQuantABA) {
				this.quantABANew = minAdmisQuantABA;
				Util.log("ABA lower bound exceeded: deficit = "+(minAdmisQuantABA - quant));
			}
			else if (quant>maxAdmisQuantABA) {
				this.quantABANew = maxAdmisQuantABA;
				Util.log("ABA upper bound exceeded: excess = "+(quant-maxAdmisQuantABA));
			}
			else this.quantABANew = quant;
		}
	}
	
	/**
	 * Initialize "quantSoluteNew" for a new time step
	 * @param flag
	 */
	public void initSoluteNew(int flag) {
		if (flag==1) this.quantABANew = this.quantABA;
	}
  
	/**
	 * Set the article water potential
	 * @param pot
	 */
	public void setWaterPot(double pot) { waterPot = pot; } 

	/**
	 * Set the article axial water flux
	 * @param axialFlux
	 */
	public void setAxialWaterFlux(double axialFlux){ axialWaterFlux = axialFlux; } 

	/**
	 * Set the article radial water flux
	 * @param radialFlux
	 */
	public void setRadialWaterFlux(double radialFlux){ radialWaterFlux = radialFlux; }
  
	/**
	 * Get the water potential [MPa] of the article
	 * @return the water potential
	 */
	public double getWaterPot(){ return waterPot; }

	/**
	 * Get the radial water flux of the article
	 * @return the radial water flux
	 */
	public double getRadialWaterFlux(){ return radialWaterFlux; }

	/**
	 * Get the axial water flux of the article
	 * @return the axial water flux
	 */
	public double getAxialWaterFlux(){ return axialWaterFlux; }
}
