
package maize;

import java.util.ArrayList;
import java.util.List;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperSymmPackMatrix;
import no.uib.cipr.matrix.Vector;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

/**
 * Class solving different functions in the model: water movement, ABA and carbon
 * 
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Xavier Draye - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Vincent LArondelle - Université catholique de Louvain - Earth and Life Institute (Belgium) *  */

public class Solver {
   
	ArrayList<Node> nodeList = new ArrayList<Node>();
	ParameterSet parms;

	/**
	 * Solve the water flux in the plant with the resolution of Doussan 98
	 * @param reseau
	 * @return true if the solve method succeed
	 */
	public boolean solveWater(Network reseau, ParameterSet p) {
	   
		parms = p;
		nodeList.clear();
      
		// Read the network structure and size A and b matrix/vector accordingly
		inscrit(reseau, -1);
		int nNode = nodeList.size(); 
		Matrix A = new UpperSymmPackMatrix(nNode);
		Vector b = new DenseVector(nNode);      

		// Iterate over the space of solutions
      	
		int i = 0;
		for (Node n : nodeList) {
			int j = n.indicePere;
           
            // ArtVirt is the first article, hence defined by its children
            if (j == -1) {
               i++;
               continue; 
            }
            
            Article a = n.article;
                        
            double surface = a.getSurface() / 1.0e4; 			// [m2]
            double dist = a.getLength() / 1.0e2;				// [m]
            double Cax = a.getWaterAxialConductance();			// [m4*s-1*Mpa-1]
            double Crad = a.getWaterRadialConductance(); 		// [m*s-1*MPa-1]
            double waterPotExo = a.getWaterPotExo(); 			// [MPa]
                                    
            double alpha = Crad * surface;
            double beta = Cax / dist;
            double gamma = waterPotExo * Crad * surface;
          
            A.set(i, i, beta+alpha);
            A.add(j, j, beta);
            A.set(j, i, -beta);
            b.set(i, gamma);
            
            i++;
         }
                   
         // Solve the equation Ax = B   
         Vector x = b.copy();
         
         try {A.solve(b, x);}
         catch (Exception e) { Util.logException(e);}
         
         // Return potential to the articles and calcul de fluxes 
         i = 0;
         for (Node n : nodeList){
        	 n.article.setWaterPot(x.get(i++));
         }    
         return true;
	}
   
	/**
	 * Solve the water flux in the plant with the resolution of Doussan 98
	 * Method using the cern.colt package
	 * @param parms
	 * @return true if the solve method succeed
	 */
	public boolean solveWater2(Network reseau) {
	   
		nodeList.clear();
	  
		// Read the network structure and size A and b matrix/vector accordingly
		inscrit(reseau, -1);
		int nNode = nodeList.size();      
		DoubleMatrix2D m = new SparseDoubleMatrix2D(nNode, nNode);
		DoubleMatrix2D v = new DenseDoubleMatrix2D(nNode, 1);    
	  
		// Iterate over the space of solutions  	
		int i = 0;
		for (Node n : nodeList) {
			int j = n.indicePere;
	       
	        // ArtVirt is the first article, hence defined by its children
			if (j == -1) {
				i++;
				continue; 
			}
	        
	        Article a = n.article;
	                    
	        double surface = a.getSurface() / 1.0e4; 				// [m2]
	        double dist = a.getLength() / 1.0e2;				// [m]
	        double Cax = 1.0 / a.getWaterAxialResistance();		// [m4*s-1*Mpa-1]
	        double Crad = 1.0 / a.getWaterRadialResistance(); 	// [m*s-1*MPa-1]
	        double waterPotExo = a.getWaterPotExo(); 			// [MPa]
	        
	        double alpha = Crad * surface;
	        double beta = Cax / dist;
	        double gamma = waterPotExo * Crad * surface;
	     
	        m.set(i, i, beta+alpha);
	        m.set(j, j, m.get(j, j) + beta);
	        m.set(j, i, -beta);            
	        m.set(i, j, -beta);            
	        v.set(i, 0, gamma);        
	        i++;
	     }
	               
	     // Solve the equation Ax = B   
	     Algebra al = new Algebra();
	     DoubleMatrix2D v2 = al.solve(m, v);
	
	     
	     // Return potential to the articles and calcul de fluxes     
	     i = 0;
	     for (Node n : nodeList){
	    	 n.article.setWaterPot(v2.get(i++, 0));
	     }    
	     return true;
	}
	
	
	
	/**
    * Solves the solute flows in the xylem
    * This is inspired by the solveWater method, adapted to solve the solute flows.
    * @author Vincent Larondelle (UCL)
    * @param reseau
    * @param parms
    * @param flag
    */
	public void solveSoluteXylem(Network reseau, ParameterSet parms, int flag) {

		nodeList.clear();
		this.parms = parms;
	      
		// Read the network structure and size A and b matrix/vector accordingly
		inscrit(reseau, -1);
		int nNode = nodeList.size(); 					// This list of node might be created within the constructor of nodes ?
		Matrix A = new UpperSymmPackMatrix(nNode+1); 	// Coefficient matrix (liters)
		Vector b = new DenseVector(nNode+1); 			// Vector of the solute quantities of each element
		Vector x = b.copy(); 							// Vector of the solutions (concentrations on the next time step)
	
		// Iterate over the space of solutions
		int i = 0; 
		for (Node n : nodeList) {
			int j = n.indicePere;
	            
			// ArtVirt is only transporter (defined by child networks)
			if (j == -1) {
				i++;
				continue;
			}
	            
			Article a = n.article;   
			double quant = a.envEndo.getQuantSolute(flag+1); 		// Loads the "new" quantity (after reaction) of the flagged solute in this article
			double waterFlow = a.getAxialWaterFlux() * 3.6e6; 		// The coefficient matrix uses the water flows [m³/s] converted in liters [l] during a defined time step.
			double volume = a.getVolume() * 10e-3; 					// Volume [cm³] converted to liters [L]
			double timeStep = Time.getTimeStep(); 					// [h] We recommend using 1 as time step for flows modeling
	            
			b.set(i,quant); 										// Vector b compiles the solute quantity in each article 
			A.add(i,i,volume); 										// On the diagonal, one find 1 + sumOutputs * timeStep
		            
			if (j == -1) {
				i++;
				continue;
			}
			
			//If the solute flows toward the father elements, the transport is an output,
			//calculated with de concentration of the element i
			if (waterFlow >= 0){									
				A.add(i, i, waterFlow * timeStep); 						// Output toward the father
				A.add(j, i, -waterFlow * timeStep); 					// Input in the father
			}
			// If the flow comes from the father, one will use the father's concentration,
			// and so add an input on the column of the father
			else if (waterFlow<0){ 
				A.add(i, j, waterFlow * timeStep); 					// Input from the father
				A.add(j, j, -waterFlow * timeStep); 				// Output from the father
			}
		}
		i++;
	
	         
		// Solve the equation A.x = b   
		try {A.solve(b, x);}
		catch (Exception e) { }
	         
		// Update the "newVariable" and sets the "flow" variable
		// flag+1 means the "newVariable" (ex : newConcABA)
		i = 0;
		for (Node n : nodeList) {
			i++;	 
			n.article.setQuantSolute((x.get(i)), flag+1);	
		}
	}
   
   
	/**
	 * Solve the carbon related functions in the plant (production, demand, allocation)
	 * @param reseau
	 * @param parms
	 */
	public void solveCarbon(Network reseau, ParameterSet p){
		nodeList.clear();
		parms = p;
		double carbonProduced = 0;
		double stemDryMass = 0;
		double leafDryMass = 0;
		double rootDryMass = 0;
		double leafMaintenanceDemand = 0;
		double stemMaintenanceDemand = 0;
		double rootMaintenanceDemand = 0;
		double leafGrowthDemand = 0;
		double stemGrowthDemand = 0;
		double rootGrowthDemand = 0;
		double waterPot = 0;
		double rootShootPartitionning = parms.rootShootPartitionning;
		ArtVirt seed = null;
		
		// Read the network structure 
		inscrit(reseau, -1);
	      	
	      	
		// Iterate over the space of solutions
		// compute the total amout of carbohydrate, based on the quantity of water transpired (Clausnitzer94)
		for (Node n : nodeList) {
			if(n.article instanceof SegLeaf){
				SegLeaf l = (SegLeaf) n.article;
				carbonProduced +=  l.getPhotosynthesisFarquar();
				leafMaintenanceDemand += l.getMaintenanceDemand();
				leafGrowthDemand += l.getGrowthDemand();
				leafDryMass += n.article.dryMass;
			}
			else if(n.article instanceof MerisLeaf){
				leafMaintenanceDemand += n.article.getMaintenanceDemand();
				leafGrowthDemand += n.article.getGrowthDemand();
			}
			else if(n.article instanceof SegRoot){
				rootMaintenanceDemand += n.article.getMaintenanceDemand();
				rootDryMass += n.article.dryMass;
				rootGrowthDemand += n.article.getGrowthDemand();

			}
			else if(n.article instanceof MerisRoot){
				rootMaintenanceDemand += n.article.getMaintenanceDemand();
				rootGrowthDemand += n.article.getGrowthDemand();

			}
			else if(n.article instanceof SegStem){
				stemMaintenanceDemand += n.article.getMaintenanceDemand();
				stemGrowthDemand += n.article.getGrowthDemand();
				stemDryMass += n.article.dryMass;
			}
			else if(n.article instanceof MerisStem){
				stemMaintenanceDemand += n.article.getMaintenanceDemand();
				stemGrowthDemand += n.article.getGrowthDemand();
			}
			else if(n.article instanceof ArtVirt){
				seed = (ArtVirt) n.article;
				double seedCarbon = seed.getSeedCarbon();
//				System.out.println("seed = "+seedCarbon);
				carbonProduced += seedCarbon;
				waterPot = seed.getWaterPot();
			}
		}	     
	    
		
		
		double totalDryMass = stemDryMass + leafDryMass + rootDryMass;
		double reserveContrib = parms.reserveC * parms.reserveSupplyRate * Time.getTimeStep();
		carbonProduced += reserveContrib;
		parms.reserveC -= reserveContrib;
		
		// Evolutive root to shoot partitionning
		if(parms.variablePartitionning) rootShootPartitionning = 0.333 - 0.444 * waterPot;
		if(rootShootPartitionning > 1) rootShootPartitionning = 1;
	      
//		System.out.println(rootShootPartitionning);
		
		// C available for the different organs
		double carbonGrowth = (carbonProduced - (leafMaintenanceDemand + rootMaintenanceDemand + stemMaintenanceDemand))/10;
		if(carbonGrowth < 0) carbonGrowth = 0;
		double carbonShoot = carbonGrowth * (1 - rootShootPartitionning);
		double carbonStem = carbonShoot - leafGrowthDemand;
		double carbonShootLeft = carbonStem - stemGrowthDemand;
//		double leafC = (1 + (carbonShoot - leafGrowthDemand) / leafGrowthDemand);
		double leafC = carbonShoot / leafGrowthDemand;
		if(leafGrowthDemand == 0) leafC = 0;
//		double stemC = (1 + (carbonStem - stemGrowthDemand) / stemGrowthDemand);
		double stemC = carbonStem / stemGrowthDemand;
		double carbonRoot = (carbonGrowth  *  rootShootPartitionning);// + (carbonShootLeft);		
		double rootC = carbonRoot / rootGrowthDemand;
		
//		System.out.println("Carbon growth = "+carbonGrowth);
//		System.out.println("Carbon growth 2 = "+(carbonRoot + carbonShoot));
//		System.out.println("Carbon root = "+carbonRoot);
//		System.out.println("Carbon shoot = "+carbonShoot);
//		System.out.println("Carbon growth = "+carbonGrowth);
//		System.out.println("LeafDemand = "+leafGrowthDemand);
//		System.out.println("LeafC = "+leafC);
//		System.out.println("rootC = "+rootC);
//		System.out.println("Partitioning = "+rootShootPartitionning);
//		System.out.println("----------");
		
		// Create an ordered list of the roots segment and meristems based on their growth demand		
		List<Article> lRootSeg = new ArrayList<Article>();
		List<Article> lRootMeris = new ArrayList<Article>();
		for (Node n : nodeList) {
			if(n.article instanceof SegRoot){
				if(n.article.getGrowthDemand() > 0){
					int i = 0;
					while(i < lRootSeg.size() && lRootSeg.get(i).getGrowthDemand() > n.article.getGrowthDemand()) i++;
					lRootSeg.add(i, n.article);
				}
			}
		}
		for (Node n : nodeList) {
			if(n.article instanceof MerisRoot){
				if(n.article.getGrowthDemand() > 0){
					int i = 0;
					while(i < lRootMeris.size() && lRootMeris.get(i).getGrowthDemand() > n.article.getGrowthDemand()) i++;
					lRootMeris.add(i, n.article);
				}
			}
		}
		
		// Allocate the C available based on the root ordered lists
		for(int i = 0 ; i < lRootMeris.size(); i++){
			if(carbonRoot < 0) carbonRoot = 0;
			if(lRootMeris.get(i).getGrowthDemand() < carbonRoot){
				lRootMeris.get(i).setGrowthEff(1);
				carbonRoot -= lRootMeris.get(i).getGrowthDemand();
			}
			else {
				lRootMeris.get(i).setGrowthEff(carbonRoot/lRootMeris.get(i).getGrowthDemand());
				carbonRoot -= lRootMeris.get(i).getGrowthDemand();
			}
		}

		for(int i = 0 ; i < lRootSeg.size(); i++){
			if(carbonRoot < 0) carbonRoot = 0;
			if(lRootSeg.get(i).getGrowthDemand() < carbonRoot){
				lRootSeg.get(i).setGrowthEff(1);
				carbonRoot -= lRootSeg.get(i).getGrowthDemand();
			}
			else {
				lRootSeg.get(i).setGrowthEff(carbonRoot/lRootSeg.get(i).getGrowthDemand());
				carbonRoot -= lRootSeg.get(i).getGrowthDemand();
			}
		}

		// Allocate the C in the leaves and the stem
		if(leafC > 1) leafC = 1;
		else if(leafC < 0) leafC = 0;
	         
		if(stemC > 1) stemC = 1;
		else if(stemC < 0) stemC = 0;
	         
		for (Node n : nodeList) {
			if(n.article instanceof SegLeaf || n.article instanceof MerisLeaf ) n.article.setGrowthEff(leafC);	        	 
			else if(n.article instanceof SegStem || n.article instanceof MerisStem ) n.article.setGrowthEff(stemC);	
		}	 

		// Set the reserve in the plant
		if(carbonRoot < 0) carbonRoot = 0;
		parms.reserveC += carbonRoot;
		parms.reserveCMax = totalDryMass / 50;
		if(parms.reserveC > parms.reserveCMax) parms.reserveC = parms.reserveCMax;
	}
   
	
	
	
	
	/**
	 * Solve the nitrogen uptake and distribution in the plant
	 * @param reseau
	 * @param p
	 */
	public void solveNitrogen(Network reseau, ParameterSet p){
		
		// Try to add an N limitation function
		double Vmaxb1 = 0;
		double alphab1 = 0;
        double kpLN = 0.2;
        double lnb0 = -5;
        double lnb1 = 18; 
		double LeafN_0 = 5;
		double LeafN_0_pot = 8; /* Needed for the potential curve for the nitrogen content --GL */
	 
		double LeafN = LeafN_0; /* Need to set it because it is used by CanA before it is computed */
		double LeafN_pot = LeafN_0_pot; /* Need to set it because it is used by CanA before it is computed  --GL */
		double LeafN_crit = LeafN_0; /* Need to set it because it is used by CanA before it is computed  --GL */

		double kLN = 0.4; 

		// Get the different article biomass
		double leafBiomass = 0;
		double rootBiomass = 0;
		double stemBiomass = 0;
		double newLeafBiomass = 0;
		double newRootBiomass = 0;
		double newStemBiomass = 0;		
		for (Node n : nodeList) {
			if(n.article instanceof SegLeaf || n.article instanceof MerisLeaf){
				leafBiomass +=  n.article.getDryMass();
				newLeafBiomass +=  n.article.getGrowthDemand();
			}
			else if(n.article instanceof SegRoot || n.article instanceof MerisRoot){
				rootBiomass +=  n.article.getDryMass();
				newRootBiomass +=  n.article.getGrowthDemand();
			}
			else if(n.article instanceof SegStem || n.article instanceof MerisStem){
				stemBiomass +=  n.article.getDryMass();
				newStemBiomass +=  n.article.getGrowthDemand();
			}
		}
		
		
//				//LeafN = LeafN_0 * exp(-kLN * TTc);				// Current nitrogen content  --GL
//
//		LeafN_crit = LeafN_0 * Math.pow(leafBiomass + stemBiomass,-kLN);		 // Critical nitrogen content --GL
//		LeafN_pot = LeafN_0_pot * Math.pow(leafBiomass + stemBiomass,-kLN);		// Potential Nitrogen content --GL
//
//		double DemandN = LeafN_0_pot * (newLeafBiomass + newStemBiomass);
//
//				// Define the N in the soil
//				
//		        /* Nitrogen fertilizer */
//		        /* Only the day in which the fertilizer was applied this is available */
//				/* When the day of the year is equal to the day the N fert was applied
//		 		/ * then there is addition of fertilizer */
//				if(doyNfert == *(pt_doy+i)){
//					Nfert = REAL(CENTCOEFS)[17] / 24.0;
//				}else{
//					Nfert = 0;
//				} 
//				minNitroPot = minNitroPot + Nfert;
//
//		 
//				NConcPot = minNitroPot * VolumePot; // Carefull with units! needs to be in micro moles / litres  // TODO
//
//				// Define N uptake capacity
//				double vmaxN1 = 0.0018;
//				double kmaxN1 = 50;
//				double vmaxN2 = 0.05;
//				double kmaxN2 = 25000;
//
//				UptakeCapacityN = (vmaxN1 * NConcPot) / (kmaxN1 + NConcPot) + (vmaxN2 * NConcPot) / (kmaxN2 + NConcPot)  ;
//				
//				// N fluxes in the roots 
//				RootFluxN = 33.6 * UptakeCapacityN * Root * rootSpecificLength; // TODO check 33.6
//
//				// Define the N supply
//				SupplyN = min(RootFluxN, minNitroPot);
//
//
//				double ratioN = min(1, (DemandN / SupplyN));
//
//				double AbsorbedN = ratioN * SupplyN;
//
//				minNitroPot = minNitroPot - AbsorbedN;
//
//				LeafN = LeafN + (AbsorbedN / (Leaf+Stem));
//
//				// Stress function
//				double LeafNS = min(1, (LeafN / LeafN_crit));
//
//
//				// Effect on the total assimilation
//
//				CanopyA = CanopyA * LeafNS;
		
		
	}
	
	
	
   /**
    * Reset the growthEff value for all the articles  
    * @param reseau
    */
	public void resetGrowthEff(Network reseau){
		nodeList.clear();
		inscrit(reseau, -1);
		for (Node n : nodeList) {n.article.setGrowthEff(1);}
	}
	        
	/**
	 * Create a node list with all the article forming the plant
	 * @param r
	 * @param pere
	 */
	public void inscrit(Network r, int pere) {
		addNode(r, pere);
		int np = nodeList.size() - 1;
		if (r.childNetworkList != null) 
			for (Network rf : r.childNetworkList) 
				inscrit(rf, np);
	}
   
	/**
	 * Add node the the node list
	 * @param r
	 * @param pere
	 */
	public void addNode(Network r, int pere) {
		nodeList.add(new Node(r.baseArticle, pere));
	}

	/**
	 * Node class used for the node list containing all the articles
	 * @author Xavier Draye - Université catholique de Louvain - Earth and Life Institute (Belgium)
	 */
	class Node {
		public Article article;
		public int indicePere;
		public Node(Article a, int i) {
			article = a;
			indicePere = i;
		}
	}
}



//import cern.colt.matrix.*;
//import cern.colt.matrix.impl.DenseDoubleMatrix2D;
//import cern.colt.matrix.impl.SparseDoubleMatrix2D;
//import cern.colt.matrix.linalg.Algebra;

//	/**
//	 * Solve the water flux in the plant with the resolution of Doussan 98
//	 * Method using the cern.colt package
//	 * @param reseau
//	 * @param parms
//	 * @return true if the solve method succeed
//	 */
//	public boolean solveWater2(Network reseau, PlanetCTParameterSet parms) {
//	   
//		nodeList.clear();
//		this.parms = parms;
//	  
//		// Read the network structure and size A and b matrix/vector accordingly
//		inscrit(reseau, -1);
//		int nNode = nodeList.size();      
//		DoubleMatrix2D m = new SparseDoubleMatrix2D(nNode, nNode);
//		DoubleMatrix2D v = new DenseDoubleMatrix2D(nNode, 1);    
//	  
//		// Iterate over the space of solutions  	
//		int i = 0;
//		for (Node n : nodeList) {
//			int j = n.indicePere;
//	       
//	        // ArtVirt is the first article, hence defined by its children
//			if (j == -1) {
//				i++;
//				continue; 
//			}
//	        
//	        Article a = n.article;
//	                    
//	        double surface = a.surface() / 1.0e4; 				// [m2]
//	        double dist = a.getLength() / 1.0e2;				// [m]
//	        double Cax = 1.0 / a.getWaterAxialResistance();		// [m4*s-1*Mpa-1]
//	        double Crad = 1.0 / a.getWaterRadialResistance(); 	// [m*s-1*MPa-1]
//	        double waterPotExo = a.getWaterPotExo(); 			// [MPa]
//	        
//	        double alpha = Crad * surface;
//	        double beta = Cax / dist;
//	        double gamma = waterPotExo * Crad * surface;
//	     
//	        m.set(i, i, beta+alpha);
//	        m.set(j, j, m.get(j, j) + beta);
//	        m.set(j, i, -beta);            
//	        m.set(i, j, -beta);            
//	        v.set(i, 0, gamma);        
//	        i++;
//	     }
//	               
//	     // Solve the equation Ax = B   
//	     Algebra al = new Algebra();
//	     DoubleMatrix2D v2 = al.solve(m, v);
//	
//	     
//	     // Return potential to the articles and calcul de fluxes     
//	     i = 0;
//	     for (Node n : nodeList){
//	    	 n.article.setWaterPot(v2.get(i++, 0));
//	     }    
//	     return true;
//	}
//	
//	/**
//	 * Solve the carbon related functions in the plant (production, demand, allocation)
//	 * @param reseau
//	 * @param parms
//	 */
//	public void solveCarbon(Network reseau, PlanetCTParameterSet parms){
//		nodeList.clear();
//		this.parms = parms;
//		double carbon = 0;
//		double stemDryMass = 0;
//		double leafDryMass = 0;
//		double rootDryMass = 0;
//		double leafMaintenanceDemand = 0;
//		double stemMaintenanceDemand = 0;
//		double rootMaintenanceDemand = 0;
//		double leafGrowthDemand = 0;
//		double stemGrowthDemand = 0;
//		double primRootGrowthDemand = 0;
//		double secRootGrowthDemand = 0;
//		ArtVirt seed = null;
//		// Read the network structure 
//		inscrit(reseau, -1);
//	      	
//		double seedContrib = 0;
//	      	
//		// Iterate over the space of solutions
//		// compute the total amout of carbohydrate, based on the quantity of water transpired (Clausnitzer94)
//		for (Node n : nodeList) {
//			if(n.article instanceof SegLeaf){
//				SegLeaf l = (SegLeaf) n.article;
//				carbon +=  l.getPhotosynthesis();//parms.wue * (Math.abs(l.getRadialWaterFlux())  * parms.valeurPas * 3600 * 1e6);
//				leafMaintenanceDemand += l.getMaintenanceDemand();
//				leafGrowthDemand += l.getGrowthDemand();
//				leafDryMass += n.article.dryMass;
//			}
//			else if(n.article instanceof MerisLeaf){
//				leafMaintenanceDemand += n.article.getMaintenanceDemand();
//				leafGrowthDemand += n.article.getGrowthDemand();
//			}
//			else if(n.article instanceof SegRoot){
//				rootMaintenanceDemand += n.article.getMaintenanceDemand();
//				if(n.article.getOrdre() == 1) primRootGrowthDemand += n.article.getGrowthDemand();
//				else secRootGrowthDemand += n.article.getGrowthDemand();
//				rootDryMass += n.article.dryMass;
//			}
//			else if(n.article instanceof MerisRoot){
//				rootMaintenanceDemand += n.article.getMaintenanceDemand();
//				primRootGrowthDemand += n.article.getGrowthDemand();
//			}
//			else if(n.article instanceof SegStem){
//				stemMaintenanceDemand += n.article.getMaintenanceDemand();
//				stemGrowthDemand += n.article.getGrowthDemand();
//				stemDryMass += n.article.dryMass;
//			}
//			else if(n.article instanceof MerisStem){
//				stemMaintenanceDemand += n.article.getMaintenanceDemand();
//				stemGrowthDemand += n.article.getGrowthDemand();
//			}
//			else if(n.article instanceof ArtVirt){
//				seed = (ArtVirt) n.article;
//				seedContrib = seed.getSeedCarbon();
//				carbon += seedContrib;
//			}
//		}	     
//	         
//		double totalDryMass = stemDryMass + leafDryMass + rootDryMass;
//		double reserveContrib = reserveC * reserveSupplyRate * Time.getTimeStep();
//		carbon += reserveContrib;
//		reserveC -= reserveContrib;
//	         
//	//	double carbonGrowth = carbon - (leafMaintenanceDemand + rootMaintenanceDemand + stemMaintenanceDemand);
//		double carbon2 = (carbon - (leafMaintenanceDemand + rootMaintenanceDemand + stemMaintenanceDemand));
//		double carbonShoot = carbon2 * (1 - parms.rootShootPartitionning);
//		double carbonStem = carbonShoot - leafGrowthDemand;
//	//	double carbonGrowth = carbon2 * 0.6;//carbonShoot - (leafGrowthDemand + stemGrowthDemand) + (carbon2 * 0.6);
//		double carbonShootLeft = carbonStem - stemGrowthDemand;
//		double carbonGrowth = carbon2 * parms.rootShootPartitionning + (carbonShootLeft);
//		double carbonPrimRoot = carbon2 * parms.rootShootPartitionning + (carbonShootLeft);
//		double carbonSecRoot = carbonPrimRoot - primRootGrowthDemand;
//		double carbonRootLeft = carbonSecRoot - secRootGrowthDemand;
//		double carbonLeft = carbonRootLeft + carbonShootLeft;
//		double shootC = (1 + (carbonShoot - (leafGrowthDemand + stemGrowthDemand)) / (leafGrowthDemand + stemGrowthDemand));
//		double leafC = (1 + (carbonShoot - leafGrowthDemand) / leafGrowthDemand);
//		double stemC = (1 + (carbonStem - stemGrowthDemand) / stemGrowthDemand);
//		double primRootC = (1 + (carbonPrimRoot - primRootGrowthDemand) / primRootGrowthDemand);
//		double secRootC = (1 + (carbonSecRoot - secRootGrowthDemand) / secRootGrowthDemand);
//	         
//	         
//	
//		/*
//		 * Create an order list of the root based on their growth demand
//		 */
//		List<Article> lRootSeg = new ArrayList<Article>();
//		List<Article> lRootMeris = new ArrayList<Article>();
//		for (Node n : nodeList) {
//			if(n.article instanceof SegRoot){
//				if(n.article.getGrowthDemand() > 0){
//					int i = 0;
//					while(i < lRootSeg.size() && lRootSeg.get(i).getGrowthDemand() > n.article.getGrowthDemand()) i++;
//					lRootSeg.add(i, n.article);
//				}
//			}
//		}
//		for (Node n : nodeList) {
//			if(n.article instanceof MerisRoot){
//				if(n.article.getGrowthDemand() > 0){
//					int i = 0;
//					while(i < lRootMeris.size() && lRootMeris.get(i).getGrowthDemand() > n.article.getGrowthDemand()) i++;
//					lRootMeris.add(i, n.article);
//				}
//			}
//		}
//		int unsat = 0;
//		for(int i = 0 ; i < lRootMeris.size(); i++){
//			if(carbonGrowth < 0) carbonGrowth = 0;
//			if(lRootMeris.get(i).getGrowthDemand() < carbonGrowth){
//				lRootMeris.get(i).setGrowthEff(1);
//				carbonGrowth -= lRootMeris.get(i).getGrowthDemand();
//			}
//			else {
//				lRootMeris.get(i).setGrowthEff(carbonGrowth/lRootMeris.get(i).getGrowthDemand());
//				carbonGrowth -= lRootMeris.get(i).getGrowthDemand();
//				unsat++;
//			}
//		}
//	
//		int unsatSeg = 0;
//		for(int i = 0 ; i < lRootSeg.size(); i++){
//			if(carbonGrowth < 0) carbonGrowth = 0;
//			if(lRootSeg.get(i).getGrowthDemand() < carbonGrowth){
//				lRootSeg.get(i).setGrowthEff(1);
//				carbonGrowth -= lRootSeg.get(i).getGrowthDemand();
//			}
//			else {
//				lRootSeg.get(i).setGrowthEff(carbonGrowth/lRootSeg.get(i).getGrowthDemand());
//				carbonGrowth -= lRootSeg.get(i).getGrowthDemand();
//				unsatSeg++;
//			}
//		}
//		carbonLeft=carbon;
//	         
//		Util.log("Unsatisfied root = "+unsatSeg+ " / unsat mersistem = "+unsat);
//	         
//	//         Util.log("ShootC = "+shootC);
//	//         if(Time.time > parms.startStress){
//	//         if(shootC > 1) {
//	//        	 SegLeaf.setGrowthEff(1);
//	//        	 SegStem.setGrowthEff(1);
//	//        	 MerisStem.setGrowthEff(1);
//	//        	 }
//	//         else if(shootC < 0) {
//	//        	 SegStem.setGrowthEff(0);
//	//        	 SegLeaf.setGrowthEff(0);
//	//        	 MerisStem.setGrowthEff(0);
//	//         }
//	//         else{
//	//        	 SegStem.setGrowthEff(shootC);
//	//        	 SegLeaf.setGrowthEff(shootC);
//	//        	 MerisStem.setGrowthEff(shootC);
//	//         }
//	////         }
//	//         
//	//         if(primRootC > 1) {SegRoot.setGrowthEffPrim(1);}
//	//         else if(primRootC < 0) SegRoot.setGrowthEffPrim(0);
//	//         else SegRoot.setGrowthEffPrim(primRootC);
//	//         
//	//         if(secRootC > 1) {
//	//        	 SegRoot.setGrowthEffSec(1);
//	//        	 MerisRoot.setGrowthEffSec(1);
//	//         }
//	//         else if(secRootC < 0){
//	//        	 SegRoot.setGrowthEffSec(0);
//	//        	 MerisRoot.setGrowthEffSec(0);
//	//         }
//	//         else {
//	//        	 SegRoot.setGrowthEffSec(secRootC);
//	//        	 MerisRoot.setGrowthEffSec(secRootC);
//	//         }
//	         
//	         /////////
//	         
//		if(primRootC > 1) primRootC = 1;
//		else if(primRootC < 0) primRootC = 0;
//	         
//		if(secRootC > 1) secRootC = 1;
//		else if(secRootC < 0) secRootC = 0;
//	         
//	//         double carbonShoot = carbonGrowth;
//	//         double shootC = (1 + (carbonShoot - (leafGrowthDemand + stemGrowthDemand)) / (leafGrowthDemand + stemGrowthDemand));
//		if(leafC > 1) leafC = 1;
//		else if(leafC < 0) leafC = 0;
//	         
//		if(stemC > 1) stemC = 1;
//		else if(stemC < 0) stemC = 0;
//	         
//		if(shootC > 1) shootC = 1;
//		else if(shootC < 0) shootC = 0;
//	         
//		for (Node n : nodeList) {
//			if(n.article instanceof SegLeaf || n.article instanceof MerisLeaf ) n.article.setGrowthEff(leafC);
//	        	 
//			else if(n.article instanceof SegStem || n.article instanceof MerisStem ) n.article.setGrowthEff(stemC);	
//	        	 
//	//        	 else if(n.article instanceof SegRoot){
//	//        		 if(n.article.getOrdre() == 1) n.article.setGrowthEff(primRootC);
//	//        		 else n.article.setGrowthEff(secRootC);
//	//        	 }
//	//        	 else if(n.article instanceof MerisRoot) n.article.setGrowthEff(secRootC);
//		}	 
//	         
//	         
//	         /////////
//	//         double carbonLeft = carbonShoot - (leafGrowthDemand + stemGrowthDemand);
//		if(carbonLeft < 0) carbonLeft = 0;
//	         
//		reserveC += carbonLeft;
//		reserveCMax = totalDryMass / 50;
//		if(reserveC > reserveCMax) reserveC = reserveCMax;
//	         
//	//         double totalGrowthDemand = leafGrowthDemand + stemGrowthDemand + primRootGrowthDemand + secRootGrowthDemand; 
//	//         double totalMaintenanceDemand = leafMaintenanceDemand+stemMaintenanceDemand+rootMaintenanceDemand;
//	//         
//	//         
//	//         Util.log("reserve = "+reserveC);
//	//         Util.log("reserveSupply = "+(reserveC * reserveSupplyRate * Time.valeurPas));
//	//         Util.log("--------------------------------");
//	//         Util.log("Carbon Total = "+carbon);
//	//         Util.log("--------------------------------");
//	//         Util.log("Leaf dry mass = " + leafDryMass);
//	//         Util.log("Stem dry mass = " + stemDryMass);
//	//         Util.log("Root dry mass = " + rootDryMass);
//	         
//	//         if (shootC < 1){
//	
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Photosynthesis contribution = "+(int)(((carbon-reserveContrib-seedContrib)/carbon)*100));
//	//	         Util.log("Reserve contribution = "+(int)((reserveContrib/carbon)*100));
//	//	         Util.log("Seed contribution = "+(int)((seedContrib/carbon)*100));
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Total dry mass = "+totalDryMass);
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Total maintenance demand = "+(int)((totalMaintenanceDemand/carbon)*100));
//	//	         Util.log("Total growth demand = "+(int)((totalGrowthDemand/carbon)*100));
//	//	         Util.log("Carbon left = "+ reserveC);
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Leaf maintenance = " + (int)((leafMaintenanceDemand/totalMaintenanceDemand)*100));
//	//	         Util.log("Stem maintenance = " + (int)((stemMaintenanceDemand/totalMaintenanceDemand)*100));
//	//	         Util.log("Root maintenance = " + (int)((rootMaintenanceDemand/totalMaintenanceDemand)*100));
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Leaf demand = " + (int)((leafGrowthDemand/totalGrowthDemand)*100));
//	//	         Util.log("Stem demand = " + (int)((stemGrowthDemand/totalGrowthDemand)*100));
//	//	         Util.log("PrimRoot demand = " + (int)((primRootGrowthDemand/totalGrowthDemand)*100));
//	//	         Util.log("SecRoot demand = " + (int)((secRootGrowthDemand/totalGrowthDemand)*100));
//	//	         Util.log("--------------------------------");
//	//	         Util.log("Shoot satisfaction = "+(int)(shootC*100));
//	//	         Util.log("PrimRoot satisfaction = "+(int)(primRootC*100));
//	//	         Util.log("SecRoot satisfaction = "+(int)(secRootC*100));
//	//	         Util.log("--------------------------------");
//	//         }
//	//         Util.log("--------------------------------");
//	//         Util.log("Carbon left = "+carbonLeft);
//	//         Util.log("LAI = "+LAI);
//	//         Util.log("Available carbon = "+carbon);
//	//         Util.log("Loss water = "+water);
//	//         Util.log("Dry Mass = "+dryMass);
//	//         Util.log("Maintenance demand = "+(shootMaintenanceDemand+rootMaintenanceDemand));
//	//         Util.log("Growth demand = "+(shootGrowthDemand+primRootGrowthDemand+secRootGrowthDemand));
//	//         Util.log("Total C = "+carbon);
//	//         Util.log("Total growth C = "+carbonGrowth);
//	//         Util.log("Growth % = "+percentage);
//	//         Util.log("Root growth % = "+percentage1);
//	//         Util.log("--------------------------------");
//	
//	}
//
///**
// * Solves the solute flows in the xylem
// * This is inspired by the solveWater method, adapted to solve the solute flows.
// * @author Vincent Larondelle (UCL)
// * @param reseau
// * @param parms
// * @param flag
// */
//	public void solveSoluteXylem(Network reseau, PlanetCTParameterSet parms, int SolveMethod, int flag) {
//
//		nodeList.clear();
//		this.parms = parms;
//	      
//		// Read the network structure and size A and b matrix/vector accordingly
//		inscrit(reseau, -1);
//		int nNode = nodeList.size(); // This list of node might be created within the constructor of nodes ?
//		Matrix A = new UpperSymmPackMatrix(nNode+1); // Coefficient matrix (liters)
//		Vector b = new DenseVector(nNode+1); // Vector of the solute quantities of each element
//		Vector delta = new DenseVector(nNode); // Method 4 vector of variation in solute
//		Vector J = new DenseVector(nNode); // Vector of quantitative flow of solute toward the father
//		Vector x = b.copy(); // Vector of the solutions (concentrations on the next time step)
//	
//	      // Iterate over the space of solutions
//	         int i = 0; // index of a node
//	         for (Node n : nodeList) {
//	            int j = n.indicePere; // index of the parent node
//	            
//	            // ArtVirt is only transporter (defined by child networks)
//	            if (j == -1) {
//	               i++;
//	               continue;
//	            }
//	            
//	            Article a = n.article;
//	            
//	            double quant = a.envEndo.getQuantSolute(flag+1); 	// Loads the "new" quantity (after reaction) of the flagged solute in this article
////	            double conc = a.getConcSolute(flag+1); 				// Loads the concentration
//	            double waterFlow = a.getAxialWaterFlux() * 3.6e6; 	// The coefficient matrix uses the water flows [m³/s] converted in liters [l] during a defined time step.
//	            double volume = a.volume()*10e-3; 					// Volume [cm³] converted to liters [L]
//	            double timeStep = Time.getTimeStep(); 					// [h] We recommend using 1 as time step for flows modeling
//	            
//	            // Construction of vector b
//	            b.set(i,quant); // Vector b compiles the solute quantity in each article 
//	 
//	            
//	            
//	            /* METHOD 1 - implicit resolution with multidirectionnal flows */
//	            if (SolveMethod == 1){
//
//		            A.add(i,i,volume); // On the diagonal, one find 1 + sumOutputs*timeStep
//		            
//		            if (j == -1) {
//			               i++;
//			               continue;
//			            }
//		            
//		            if (waterFlow >= 0){	// If the solute flows toward the father elements, the transport is an output,
//		            						//calculated with de concentration of the element i
//		            	A.add(i,i,waterFlow*timeStep); // Output toward the father
//		            	A.add(j,i,-waterFlow*timeStep); // Input in the father
//		            }
//		            else if (waterFlow<0){ 	// If the flow comes from the father, one will use the fater's concentration,
//		            						//and so add an input on the column of the fath er
//		            	A.add(i,j,waterFlow*timeStep); // Input from the father
//		            	A.add(j,j,-waterFlow*timeStep); // Output from the father
//		            	// the signs have been changed because the waterFlow is negative
//					}
//	            }
//	            
//	            /* METHOD 2, only when the direction of the flow is toward the father element */
//	            else if (SolveMethod == 2){
////		            A.add(i+1,i+1,1+waterFlow*timeStep); // Output toward the father
////	            	A.add(j+1,i+1,-waterFlow*timeStep); // Input in the father
//	            	Util.log("SolveMethod " + SolveMethod + "has not been corrected yet. Please use SolveMethod 1");
//	            }
//
//	            /* METHOD 3 - Explicit resolution with multidirectionnal flows*/	            
//	            else if (SolveMethod ==3){
////		            A.set(i+1,i+1,volume);
////		            if (waterFlow >= 0){
////		            	A.add(i+1,i+1,-waterFlow*timeStep); // output
////		            	A.add(j+1,i+1,waterFlow*timeStep); // input
////		            }
////		            else if (waterFlow <0) {
////		            	A.add(i+1,j+1,-waterFlow*timeStep); //input
////		            	A.add(j+1,j+1,waterFlow*timeStep); // output
////		            }
//		            Util.log("SolveMethod " + SolveMethod + "has not been corrected yet. Please use SolveMethod 1");
//	            }
//	            
//		        /* METHOD 4 - Iterative resolution with multidirectionnal flows
//		        *  Calculates the "quant" flow for this element */
//	            else if (SolveMethod ==4){
//
//	            	
////		            if (waterFlow>=0){
////		            	double soluteFlow = conc*waterFlow*timeStep;
////		            	delta.add(i+1,-soluteFlow); // output from article a
////		            	delta.add(j+1,soluteFlow); // input in father
////		            	J.add(i+1,soluteFlow); // from the son to the father
////		            }
////		            else if (waterFlow<0){
////		            	double soluteFlow = -a.getArtPere().getConcSolute(flag+1)*waterFlow*timeStep;
////		            	delta.add(i+1,soluteFlow); // input in article a
////		            	delta.add(j+1,-soluteFlow); // output from father
////		            	J.add(i+1,-soluteFlow); // from the father to the son
////		            }
//		            Util.log("SolveMethod " + SolveMethod + "has not been corrected yet. Please use SolveMethod 1");
//	            }
//	            else Util.log("An unexisting solute flow solving method was choosen. Please use SolveMethod 1");
//
//	            i++;
//	         }
//	         
//
//	         if (SolveMethod == 1 || SolveMethod == 2 || SolveMethod ==3){
//		         // Solve the equation A.x = b   
//		         try {A.solve(b, x);}
//		         catch (Exception e) { }
//	         }
//
//	         // Update the "newVariable" and sets the "flow" variable
//	         i = 0;
//	         for (Node n : nodeList) {
//	        	 i++;
//	        	 
//	        	 if (SolveMethod == 4){
//			         /* Sequel METHOD 4
//			         * Update the soluteNewQuantity
//			         * Export the values contained in J to flowSolute */
//	        		 
//		        	 x.set(i,b.get(i)+delta.get(i));
//		        	 n.article.setQuantFlowSolute(J.get(i),flag); // Reminds how much solute has flown from an article in direction of his father
//	        	 }
//
//	        	 n.article.setConcSolute((x.get(i)), flag+1);	// flag+1 means the "newVariable" (ex : newConcABA)
//	         }
//	   }

