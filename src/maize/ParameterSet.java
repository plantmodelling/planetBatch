package maize;

import org.w3c.dom.*;
import javax.xml.parsers.*;
import java.io.*;

/**
 * This class contain most of the parameters of the simulation and send it to the CrossTalk interface
 * Parameters can therefore be change during the simulation
 * 
 * @author Xavier Draye - Université catholique de Louvain - Earth and Life Institute (Belgium)
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public class ParameterSet {


	// Simulation parameters
	float simulationDuration;  				// Total simulation time [hour]
	float timeStep; 						// Time value of a time step [hour]
	String outputFileName; 					// Name of the output file.  
	
	// Environment parameters
	float baseTemp;							// Base temperature for maize growth and development [C]
	float airTemp = 25; 					// Air growing temperature [C]
	float soilTemp = 25;					// Soil temperature	[C]	
	float airWaterPot = -95;				// Air water potential [MPa]	
			
	// Leaf parameters
	float leafGrowthRate = 0.564f / 24;		// Leaf growth rate	[cm/Ch] (Drouet & Pagès 2003)
	float leafApparitionRate = 1;			// Leaf appariton rate modification factor [-]
			
	// Stem parameters
	float stemInitDiam = 0.4f;				// Stem initial diameter [cm]
	float stemGrowthRate = 0.41f;			// Stem growth rate		
	
	// Root parameters
	int branchingLevel;						// Level of branching of the root system (0 = only axes)
	int numberOfSeminals;					// Number of seminal roots
	int rootID;								// Identifier of the roots
	float rootTortuosity;					// Adding random tortuosity to the root (e.g. mimicking stones in the soil)
	boolean diamGrowthRelation;				// Flag for the relation between root diameter and root growth rate
	boolean diamStochasticity;				// Flag for the stochasticity in the root diameter of newly created roots
	boolean adventives;						// Flag for the creation of adventive roots
	boolean semAngStochasticity;			// Flag for stochasticity in the seminal intertion angle
	
	float primaryMaxLength;					// Maximal length of the primary roots [cm]
	float secondaryMaxLength;				// Maximal length of the secondary roots [cm]
	float tertiaryMaxLength;				// Maximal length of the tertiary roots [cm]
	
	float primaryLAUZ;						// Lenght of the unbranched zone on the primary roots [cm]
	float secondaryLAUZ;					// Length of the unbranched zone of the secondary roots [cm]
	
	float primaryInterBranch;				// Interbranch distance on the primary roots [cm]
	float secondaryInterBranch;				// Interbranch distance on the secondary roots [cm]
	
	float secondaryInsertAngle;				// Insertion angle of secondary roots [rad]
	float tertiaryInsertAngle;				// Insertion angle of tertiary roots [rad]
	
	float principalRootInitDiam;			// Principal root initial diameter [cm]
	float seminalRootInitDiam;				// Seminal root initial diameter [cm]
	float crownRootInitDiam;				// Crown root initial diameter [cm]
	
	float principalRootGravitropism;		// Principal root gravitropism [-]
	float seminalRootGravitropism;			// Seminal root gravitropism [-]
	float secondaryRootGravitropism;		// Secondary root gravitropism [-]
	float tertiaryRootGravitropism;			// Tertiary root gravitropism [-]
	
	// Water related parameters
	boolean resolveWater;					// Flag for the water resolution algorithm		
	boolean allowHydraulicLift;				// Flag for allow hydrolic lift in the soil-plant continuum		
	boolean aqp, cavitation, stomata;		// Flags for the water resolution regulation	
	float stomataMod, aqpMod, cavitationMod;// Modifiers for the water uptake regulation
	float radMod, axMod;					// Radial and axial modifiers for the root conductivities
	int resolveIteration;					// number of iteration for the water resolution between two simulation time step
	int resolveWaterType;					// Define which algorithm will be used to solve the water [1 = cipr package, 2 = cern package]
	int stressType;							// Type of water stress applyied to the plant
	int startStress;						// Time to start the stress
	float permanentWiltingPoint = -1.5f;	// Permanent Wilting Point 
	int reWateringFrequence = 50;			// Frequence of soil rewatering
	int splitRootType = 0;					// Type of split root scenario to apply to the plant
	int splitRootLimit = 0;					// Limit of the split root scenario to apply to the plant
	
	// ABA related parameters
	static final double molMassABA = 264.32;	// MM ABA [g/mol]
	boolean ABASolveMethod;						// Flag for the ABA resolution algorithm
	boolean ABAProduction;						// Flag for the ABA production algorithm
	
	// Carbon related parameters
	boolean resolveCarbon;													// Flag for the carbon resolution algorithm	
	boolean variablePartitionning;											// Flag for a carbon partitionning based on the plant water status						
	float[] wueX = {0, 0.2f, 0.4f, 0.6f, 0.8f, 1};							// X vector for the WUE-gs relationship (Allen et al 2011)
	float[] wueY = {6e-3f, 8.0e-3f, 9.2e-3f, 9.0e-3f, 7.5e-3f, 6.0e-3f};	// Y vector for the WUE-gs relationship (Allen et al 2011)
	float rootShootPartitionning; 											// Root to shoot partitioning ration for C
	double reserveC = 0;												// Carbon reserve in the plant
	double reserveCMax = 0.1f;										// Maximal carbon reserve allowed in the plant
	final float reserveSupplyRate = 0.05f;							// Carbon reserve supply rate
	float maxPAR;													// Max daily PAR
	float maxTemperature;													// Max daily temeprature [°C]
	float VMax;									// maximal rate of carboxylation
	float JMax;									// maximal rate of electron transport 
	float aerenchyma;									// proportion of aerenchyma in the roots 
	
	// Export parameters
	int exportType;								// Export type: 0 = no export, 1 = SQL, 2 = CSV
	boolean exportAtNoon;						// Flag for export only at 12:00
	boolean exportArchitecture;					// Flag for architecture export
	boolean exportWater;						// Flag for water related export
	boolean exportABA;							// Flag for ABA related export
	boolean exportCarbon;						// Flag for carbon related export	
	boolean exportResEq;						// Flag for export of the equivalent resistances
	int startExport;							// Time to start the export
	
	// Lookup table parameters
	float minPot = 1.5f, maxPot = 0;			// Max / Min value for the potential LUT
	float minPotExo = -5f, maxPotExo = 0;			// Max / Min value for the soil potential LUT
	float maxAFlux = 1.0e-9f, minAFlux = 0; 	// Max / Min value for the axial flux LUT
	float maxRFlux = 1.0e-9f, minRFlux = 0;		// Max / Min value for the radial flux LUT
	float maxABA = 1.0e3f, minABA = 0;			// Max / Min value for the ABA LUT
	float maxRRes = 2.5e7f, minRRes = 0.2e7f;	// Max / Min value for the radial resistance LUT
	float maxARes = 6.0e9f, minARes = 0;		// Max / Min value for the axial resistance LUT
	
	
	private static Document xml;
	
	
	
   /**
    * Load the XML file specified in the args[]
    */
	public ParameterSet(String file){
		try{
	        File inputFile = new File(file);
			DocumentBuilderFactory factory =
			DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			
	         xml = builder.parse(inputFile);
	         xml.getDocumentElement().normalize();
	         
		}catch(Exception e){
			System.out.println("XML document could not be loaded: "+e);
		}
		
	}
     
  
	/** @see simulator.model.ParameterSet#updateParameterValues() */
	public boolean updateParameterValues() {
		
		Util.log("Updating values of " + getClass().getName());
 
		XPathParser xpp = new XPathParser(xml.getDocumentElement());

		// Simulation	
		
		baseTemp = xpp.getFloat("Options/Simulation/BaseTemperature");
		timeStep = xpp.getInteger("Options/Simulation/TimeStep");
		simulationDuration = xpp.getInteger("Options/Simulation/TotalTime");
		outputFileName = xpp.getString("Options/Simulation/outputFileName");
		
		// Export
		exportType = xpp.getInteger("Options/Export/ExportType");
		exportAtNoon = xpp.getInteger("Options/Export/ExportAtNoon") == 1;
		exportArchitecture = xpp.getInteger("Options/Export/ExportArchitecture") == 1;
		exportWater = xpp.getInteger("Options/Export/ExportWater") == 1;
		exportABA = xpp.getInteger("Options/Export/ExportABA") == 1;
		exportCarbon = xpp.getInteger("Options/Export/ExportCarbon") == 1;
		exportResEq = xpp.getInteger("Options/Export/ExportResEq") == 1;
		startExport = xpp.getInteger("Options/Export/StartExport");
		
		// Water
		resolveWater = xpp.getInteger("Options/Water/ResolveWater") == 1;
		splitRootType = xpp.getInteger("Options/Water/SplitRoot");
		splitRootLimit = xpp.getInteger("Options/Water/SplitRootLimit");
		allowHydraulicLift = xpp.getInteger("Options/Water/AllowHydrolicLift") == 1;
		resolveWaterType = xpp.getInteger("Options/Water/ResolveWaterType");
		resolveIteration = xpp.getInteger("Options/Water/ResolveIteration");
		reWateringFrequence = xpp.getInteger("Options/Water/ReWateringFrequence");
		aqp = xpp.getInteger("Options/Water/Aquaporines") == 1;
		cavitation = xpp.getInteger("Options/Water/Cavitation") == 1;
		stomata = xpp.getInteger("Options/Water/Stomata") == 1;
		stressType = xpp.getInteger("Options/Water/StressType");
		startStress = xpp.getInteger("Options/Water/StartStress");
		radMod = xpp.getFloat("Options/Water/RadialModifier");
		axMod = xpp.getFloat("Options/Water/AxialModifier");
		stomataMod = xpp.getFloat("Options/Water/Regulation/StomataModifier");
		aqpMod = xpp.getFloat("Options/Water/Regulation/AQPModifier");
		cavitationMod = xpp.getFloat("Options/Water/Regulation/CavitationModifier");

		// ABA
		ABASolveMethod = xpp.getInteger("Options/ABA/ABATransport") == 1;
		ABAProduction = xpp.getInteger("Options/ABA/ABAProduction") == 1;

		// Carbon
		resolveCarbon = xpp.getInteger("Options/Carbon/ResolveCarbon") == 1;
		variablePartitionning = xpp.getInteger("Options/Carbon/VariablePartitionning") == 1;
		rootShootPartitionning = xpp.getFloat("Options/Carbon/RootToShootPartition");
		maxPAR = xpp.getFloat("Options/Carbon/maxPAR");
		maxTemperature = xpp.getFloat("Options/Carbon/maxTemperature");
		VMax = xpp.getFloat("Options/Carbon/VMax");
		JMax = xpp.getFloat("Options/Carbon/JMax");
		aerenchyma =  xpp.getFloat("Options/Carbon/Aerenchyma");
		
		// Roots
		rootTortuosity  = xpp.getFloat("Root/General/RootTortuosity");
		diamStochasticity  = xpp.getInteger("Root/General/DiameterStochasticity") == 1;
		branchingLevel = xpp.getInteger("Root/General/RamificationLevel");
		numberOfSeminals = xpp.getInteger("Root/General/NumberOfSeminals");
		adventives = xpp.getInteger("Root/General/AdventiveRoots") == 1; 
		semAngStochasticity = xpp.getInteger("Root/General/SeminalAngleStochasticity") == 1;
		
		primaryMaxLength = xpp.getFloat("Root/MaxLength/PrimaryMaxLength");		
		secondaryMaxLength = xpp.getFloat("Root/MaxLength/SecondaryMaxLength");
		tertiaryMaxLength = xpp.getFloat("Root/MaxLength/TertiaryMaxLength");

		primaryLAUZ = xpp.getFloat("Root/LAUZ/PrincipalLAUZ");
		secondaryLAUZ = xpp.getFloat("Root/LAUZ/SecondaryLAUZ");
		
		primaryInterBranch = xpp.getFloat("Root/InterBranchDistances/PrimaryInterBranchDistance");		
		secondaryInterBranch = xpp.getFloat("Root/InterBranchDistances/SecondaryInterBranchDistance");	
		
		secondaryInsertAngle = xpp.getFloat("Root/InsertAngles/SecondaryInsertAngle");	
		tertiaryInsertAngle = xpp.getFloat("Root/InsertAngles/TertiaryInsertAngle");	
		
		principalRootInitDiam = xpp.getFloat("Root/Diameters/PrincipalInitDiameter");	
		seminalRootInitDiam = xpp.getFloat("Root/Diameters/SeminalInitDiameter");	
		crownRootInitDiam = xpp.getFloat("Root/Diameters/CrownInitDiameter");	
		
		principalRootGravitropism = xpp.getFloat("Root/Gravitropisms/PrincipalRootGravitropism");	
		seminalRootGravitropism = xpp.getFloat("Root/Gravitropisms/SeminalRootGravitropism");	
		secondaryRootGravitropism = xpp.getFloat("Root/Gravitropisms/SecondaryRootGravitropism");	
		tertiaryRootGravitropism = xpp.getFloat("Root/Gravitropisms/TertiaryRootGravitropism");	

		
		// Shoot
		leafGrowthRate = xpp.getFloat("Shoot/Leaf/LeafGrowthRate") / 24;
		leafApparitionRate = xpp.getFloat("Shoot/Leaf/LeafApparitionRate");
		return true;
   }

	/** @see simulator.model.ParameterSet#validate() */
	public boolean validate() {
		Util.log("Validation of model " + getClass().getName());
		return true;
	}
}
