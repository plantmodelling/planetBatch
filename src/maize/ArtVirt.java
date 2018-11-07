package maize;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Virtual article, as a basis to the network (position, flux)
 * represented in the form of a cylinder with a fixed length and a diameter
 * The virtual article is only transporter
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium) 
 */

public class ArtVirt extends Article {
	
	public final double diameter=1.0;
	public static final float seedSupplyRate = 0.3f; 		// Seed supply rate of carbon. From Drouet and Pagès 2003
	public static final float initialSeedDryMass = 0.5f;	// Seed initial dry mass. From Drouet and Pagès 2003
	public float seedDryMass = initialSeedDryMass;			// Seed current dry mass.
  
	/**
	 * Constructor
	 * @param node
	 * @param dirAx
	 * @param dirNorm
	 * @param envEndo
	 * @param parms
	 */
	public ArtVirt(Point3d node, Vector3d dirAx, Vector3d dirNorm, EndogenousEnvironment envEndo, ParameterSet parms) {
		super(node, dirAx, dirNorm, envEndo, parms);
	}
	
	/**
	 * Develop the article. In this case, the article is not growing
	 */
	public void develop() {
		updateThermalAge();
	}

	////////////////////////////////////  GET FUNCTIONS  /////////////////////////////////////////
	
	/**
	 * Get the volume of the article, considered as a cylinder 
	 * @category architecture
	 * @return the volume [cm3]
	 */
	public double getVolume() { return Util.getSphereVolume(diameter);}
  
	/**
	 * Get the surface of the article, considered as a cylinder 
	 * @category architecture
	 * @return the surface [cm2]
	 */
	public double getSurface(){ return Util.getSphereSurface(diameter); }
	
	/**
	 * Get the length of the article. As this article is round, return the diameter
	 * @category architecture
	 * @return the length
	 */
	public double getLength() { return diameter; }

	/**
	 * Get the width of the article. As this article is round, return the diameter
	 * @category architecture
	 * @return the width
	 */
	public double getWidth() { return diameter;}

	/**
	 * Get the thickness of the article. As this article is round, return the diameter
	 * @category architecture
	 * @return the thickness
	 */
	public double getThickness() { return diameter;}

	/**
	 * Get the ABA consumption of the article
	 * Zero in the virtual article
	 * @author Vincent Larondelle
	 * @category solute
	 * @return 0  
	 */
	public double getABAConsumption() { return 0; }

	/**
	 * Get the ABA consumption of the article
	 * Zero in the virtual article
	 * @author Vincent Larondelle
	 * @category solute
	 * @return 0  
	 */
	public double getABAProduction() {	return 0;}
 
	/**
	 * Return the axial resistance of the article
	 * The axial resistance of the virtual article is the average of the resistances of its leaf and root children
	 * @category water
	 * @return axial resistance
	 */
	public double getWaterAxialResistance() {
		double khRoot = 0;
		int nRoot = 0;
	  	double khStem = 0;
	  	int nLeaf = 0;
	  	double kh = 0;
		if (network.childNetworkList != null) {
			for (Network resFils : network.childNetworkList) {
				if(resFils.baseArticle instanceof SegRoot){
					khRoot += resFils.baseArticle.getWaterAxialResistance();
					nRoot ++;
				}
				else if(resFils.baseArticle instanceof SegStem){
					khStem += resFils.baseArticle.getWaterAxialResistance();
					nLeaf++;
				}
			}
		}
		kh = ((khRoot/nRoot)+(khStem/nLeaf))/2;
		
	    return kh;
	}
  
	/**
	 * Get the radial resistance of the article
	 * Radial conductance of the artvir is 0
	 * @category water
	 * @return the radial resistance
	 */
    public double getWaterRadialConductance(){ return Util.epsilonDouble; }
    
	/**
	 * Get the axial conductanc of the article. Zero for the virtual article
	 * @category water
	 * @return the axial conductance
	 */
	public double getWaterAxialConductance() { return 1 / getWaterAxialResistance(); }

    
    /**
     * Get the equivalent resistance of the article
     * @category water
     * @return the equivalent resistance
     */
	public double getEquivalentConductance() {
		double s = 0;
		if (network.childNetworkList != null) {
			for (Network resFils : network.childNetworkList) {
				s += (resFils.baseArticle.getLength() / 1.0e2 / (1 / getWaterAxialResistance())) + resFils.baseArticle.KhEqu;
			}
		}
		KhEqu = s;
		return s;	
	}
	
	/**
	 * Get the growth demand of the article
	 * As the virtual article does not grow, return 0
	 * @category carbon
	 * @return 0
	 */
	public double getGrowthDemand(){ return 0; }
	
	/**
	 * Get the quantity of available carbon [g] in the seed at this time step
	 * Update the carbon stock
	 * From Drouet and Pagès 2003
	 * @category carbon
	 * @return the carbon available for growth
	 */
	public double getSeedCarbon(){
		float c = seedDryMass * seedSupplyRate * (float) Time.getTimeStep();
		updateSeedStock(-c);
		return c;
	}
	
	/**
	 * Update the stock of carbon in the seed
	 * @category carbon
	 * @param c the C quantity uptaken
	 */
	public void updateSeedStock(float c){
		seedDryMass = seedDryMass + c;
		if(seedDryMass < 1e-5) seedDryMass = 0;
	}

	/**
	 * Get the growth. Zero for the virtual article
	 * @category architecture
	 * @return 0
	 */
	public float getGrowth() { return 0; }

	/**
	 * Get the distance from the base article. Zero for the virtual article
	 * @category architecture
	 * @return 0
	 */
	public float getDistFromBase() { return 0; }

	/**
	 * Get the article ID. Zero for the virtual article
	 * @category architecture
	 * @return 0
	 */
	public int getArticleID() { return 0;}

	/**
	 * Get aquaporin effect. Zero for the virtual article
	 * @category water
	 * @return 0
	 */	
	public double getAQPCondEffect() { return 0; }

	/**
	 * Get the reduction in radial conductivity. Zero for the virtual article
	 * @category architecture
	 * @return 0
	*/	 
	public float getLrPercent() { return 0; }

	/**
	 * Get the photosynthesis of the article. Zero for the virtual article
	 * @category architecture
	 * @return 0
	 */
	public double getPhotosynthesis() { return 0; }

	/**
	 * Get the maintenance demand. Zero for the virtual article
	 * @category carbon
	 * @return 0
	 */
	public double getMaintenanceDemand() { return 0;}

}
