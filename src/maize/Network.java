package maize;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point3d;

/** 
 * Represents a network of articles, that is to say a 
 * branched system (tree graph) in which the items are connected.
 * This class manages the connection between aspects of the 
 * articles, and courses (acropetal or basipetal).
 * The system aerial and root system are typically networks.
 * We can possibly consider other descendants of the network, 
 * like the seed, the metamer, or any other composed structure.
 * 
 * @author Loic Pagès, INRA Avignon
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent Larondelle - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public class Network {
	
	public Article baseArticle;        				// base article of the network
	public Network parentNetwork;         			// parent network, on which this one is attached
	public List<Network> childNetworkList;     		// list of the networks attached to this one
	static Article closestArticle = null;			// the closest article from a given point
	static double distToClosestArticle = 10000;		// the distance to the closest article
	public static int childNumber = 0;				// the number of children of the netword
  
	/**
	 * Constructor of an empty network
	 */
	public Network() {	
		baseArticle=null;
		parentNetwork=null;
		childNetworkList=null;
	}

	/**
	 * Constructor of new network, without children
	 * @param artBase
	 */
	public Network(Article artBase) {
		this.baseArticle=artBase;
		this.parentNetwork=null;
		childNetworkList=null;
		artBase.setNetwork(this);
	}  
   
	/**
	 * Get the base article of the network
	 * @return the base article
	 */
	public Article getBaseArticle() { return baseArticle; }

	/**
	 * Get the father network of the network
	 * @return the father network
	 */
	public Network getParentNetwork() { return parentNetwork; } 

	/**
	 * Get the list of the children networks
	 * @return the list of children networks
	 */
	public List<Network> getChildNetworkList() { return childNetworkList; } 


	/**
	 * Add a child network in a proximal position.
	 * The new network is positionned in the list of children based on its relative position
	 * @param network to add
	 */
	public void addProximalChildNetwork(Network network) {  

		if (childNetworkList==null) childNetworkList=new ArrayList<Network>();
		// The new network is positionned in the list of children based on its relative position    
		int i=childNetworkList.size();
		childNetworkList.add(network);
		while (i>0){
			childNetworkList.add(i+1, childNetworkList.get(i));
			i--;
		}
		childNetworkList.add(0,network); 
	}
 
	/**
	 * Add a child network in a distal position.
	 * The new network is positionned in the list of children based on its relative position
	 * @param network to add
	 */
	public void addDistalChildNetwork(Network network) {  

		if (childNetworkList==null) childNetworkList=new ArrayList<Network>();
		int i=childNetworkList.size();
		network.parentNetwork=this; // the added network take this one as father
		while ((i>0) && (childNetworkList.get(i-1).baseArticle.relativePosition > network.baseArticle.relativePosition)) i--;
		childNetworkList.add(i,network); 
	}
  
	/**
	 * Remove the specified network from the list of children
	 * @param network to remove
	 */
	public void removeChildNetwork(Network network) {  
		if (childNetworkList==null) return;
		else {
			if (!childNetworkList.remove(network)) Util.log("Network "+network.toString()+" was note removed from the list");
			if (childNetworkList.size()==0) childNetworkList=null;
		}
	}
  
	/**
	 * Remove all the network from the children list
	 */
	public void removeAllChildNetwork() {  
		if (childNetworkList==null) return;
		else {
			childNetworkList.clear();
			childNetworkList=null; 
		}   
	}
  
	/**
	 * Inserts a new network, by inserting it between the latter and his father
	 * that is to say in the position distal to the latter
	 * @param network
	 */
	public void addDistalNetwork(Network network) {  
		if(network != null){
			if(this.childNetworkList != null){
				Network n = this.childNetworkList.get(0);
				this.removeChildNetwork(n);
				this.addDistalChildNetwork(network);
				network.addDistalChildNetwork(n);
			}
			else this.addDistalChildNetwork(network);
		}
	}
  
	/**
	 * Inserts a new network, by inserting it between the latter and his father
	 * that is to say in the position proximal to the latter
	 * @param network
	 */
	public void addProximalNetwork(Network network) {  
		if (parentNetwork!=null) {
			parentNetwork.removeChildNetwork(this);        
			parentNetwork.addDistalChildNetwork(network); 
		}
		network.addDistalChildNetwork(this);
	}

	/**
	 * Positioning articles on the network, acropetally,
	 * Assumes the position of the base article is known,
	 *  from which all other positions are calculated
	 */
	public void position() {
		if (childNetworkList!=null) {  
			for (Network resFils : childNetworkList) { 
				resFils.baseArticle.setPosition(baseArticle); 
				resFils.position();  
			}
		}
	}

	/**
	 * Develop articles in the network, acropetally
	 */
	public void develop() {
	  
		baseArticle.develop(); 
		if (childNetworkList!=null) {  
			for (int i=0; i<childNetworkList.size() ; i++) {  
				Network resFils=childNetworkList.get(i);
				resFils.develop(); 
			}
		}
	}

	/**
	 * Update the value of the given solute in the network
	 * @param coeff
	 * @param flag [1 = ABA]
	 */
	public void setConcAllSolute(double conc, int flag) {
		if (childNetworkList!=null) { 
			for (Network resFils : childNetworkList) resFils.setConcAllSolute(conc, flag);
		}
		baseArticle.setQuantSolute(conc,flag);
	} 
  
	/**
	 * Update the value of the given solute in the network
	 * @param coeff
	 * @param flag [1 = ABA]
	 */
	public void setQuantSolute(double quant, int flag) {
		if (childNetworkList!=null) { 
			for (Network resFils : childNetworkList) resFils.setQuantSolute(quant, flag);
		}
		baseArticle.envEndo.setQuantSolute(quant,flag);
	} 

	/**
	 * Get the total mass of a given solute over the network
	 * @param flag [1 = ABA]
	 * @return the total mass
	 */
	public double getSoluteMass(int flag) {
		double masse=baseArticle.getMassSolute(flag);
		if (childNetworkList!=null) { 
			for (Network resFils : childNetworkList) { 
				masse+=resFils.getSoluteMass(flag);
			}
		}
		return masse;
	}

	
	/**
	 * Get the equivalent conductance of the network, basipetally
	 */
	public void getEquivalentConductance (){
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.getEquivalentConductance();
			}
		}
		baseArticle.getEquivalentConductance();
	}

	
	/**
	 * Set the radial and axial water flux in every article
	 */
	public void setWaterFlux() {
		
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.setWaterFlux();
			}
		}
		baseArticle.setRadialWaterFlux();
		baseArticle.setAxialWaterFlux();
	}
	
	/**
	 * Update the exogenous water content for every article 
	 * @param t
	 */
	public void updateWaterContentExo(double t) {
		
		baseArticle.updateWaterContExo(baseArticle.getRadialWaterFlux() * t * 1e6);
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.updateWaterContentExo(t);
			}
		}		
	}

	/**
	 * Runs soluteReaction for the whole network
	 * @author Vincent Larondelle
	 * @param flag [1 = ABA, 2 = ABA new]
	 */
	public void soluteReaction(int flag){
		baseArticle.soluteReaction(flag);
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.soluteReaction(flag);
			}
		}
	}
		
	/**
	 * Updats the quantity of a solute to it's new value
	 * @author Vincent Larondelle
	 * @param flag [1 = ABA, 2 = ABA new]
	 */
	public void updateSoluteValuesToNew(int flag) {
		baseArticle.setDeltaQuantSolute(baseArticle.envEndo.getQuantSolute(flag), baseArticle.envEndo.getQuantSolute(flag+1),flag);
		baseArticle.setDeltaConcSolute(baseArticle.getConcSolute(flag), baseArticle.getConcSolute(flag+1),flag);
		baseArticle.envEndo.setQuantSolute(baseArticle.envEndo.getQuantSolute(flag), flag); 
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.updateSoluteValuesToNew(flag);
			}
		}
	}

	/**
	 * Reset the closest article search
	 */
	public void resetArticleSearch(){
		distToClosestArticle = 10000;
		if(closestArticle != null) closestArticle.setSelected(false);
		closestArticle = null;
	}
	
	/**
	 * Get the closest article of a defined point
	 * @param p the point
	 * @return the article
	 */
	public Article getClosestArticle(Point3d p){
		if(p.distance(baseArticle.node) < distToClosestArticle){
			distToClosestArticle = p.distance(baseArticle.node);
			if(closestArticle != null) closestArticle.setSelected(false);
			closestArticle = baseArticle;
		}
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.getClosestArticle(p);
			}
		}
		closestArticle.setSelected(true);
		return closestArticle;
	}
  
	/**
	 * Get the number of children of the network
	 * @return the number of children
	 */
	public int getNumberOfChildren(){
		if(baseArticle.isSelected()) childNumber = 0;
		if (childNetworkList != null) {
			for (Network resFils : childNetworkList) {
				resFils.getNumberOfChildren();
			}
		}
		if(childNetworkList != null) childNumber += childNetworkList.size();
		return childNumber;
	}
}

