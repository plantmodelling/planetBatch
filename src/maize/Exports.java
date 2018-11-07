package maize;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.sql.Connection;
import java.sql.Statement;

import maize.Simulator.Node;

/** 
 * This class contains the exort functions to SQL database or CSV file
 * CSV export is much faster.
 * 
 * @author Guillaume Lobet - Université catholique de Louvain - Earth and Life Institute (Belgium)
 */

public class Exports {
	
	// SQL export
	static String database="planet";
	static String driver = "com.mysql.jdbc.Driver";
	static String url = "jdbc:mysql://localhost/"+database;
	static String usr = "root";
	static String psw = "";
	static String dbDate;
	static String tableND;
	static String tableParam = "simulation_parameters";
	static String tableRE;
	static String tableEnvi;

	// RSML export
	static PrintWriter rsml;	

	// CSV export
	static PrintWriter planet;
	static PrintWriter soil;
	static PrintWriter resEqu;

	static int counter = 0;
	static float[] rld = new float[70];
	static String date;
	
	/**
	 * Create RSML file which will contain the data from the simulation
	 * Create the headers.
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resequ : equivalent resistance
	 */
	public static void createRSML(){
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
		Date d = new Date();
		date = dateFormat.format(d);
		  
		try{
			rsml = new PrintWriter(new FileWriter("/Users/guillaumelobet/Desktop/"+date+"_planet.rsml"));
			rsml.println("<?xml version='1.0' encoding='UTF-8'?>\n");
			rsml.println("<rsml xmlns:po='http://www.plantontology.org/xml-dtd/po.dtd'>\n");
			rsml.println("	<metadata>\n");
			rsml.println("		<version>1</version>\n");
			rsml.println("		<unit>inch</unit>\n");
			rsml.println("		<resolution>300.0</resolution>\n");
			rsml.println("		<last-modified>today</last-modified>\n");
			rsml.println("		<software>smartroot</software>\n");
			rsml.println("		<user>globet</user>\n");
			rsml.println("		<file-key>myimage</file-key>\n");
			rsml.println("					<property-definitions>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>diameter</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>cm</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>length</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>cm</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>angle</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>degree</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>insertion</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>cm</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>lauz</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>cm</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("		    	<label>lbuz</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("				<unit>cm</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("			<property-definition>\n");
			rsml.println("				<label>node-orientation</label>\n");
			rsml.println("		        <type>float</type>\n");
			rsml.println("		        <unit>radian</unit>\n");
			rsml.println("			</property-definition>\n");
			rsml.println("		</property-definitions>\n");
			rsml.println("	</metadata>\n");
			rsml.println("	<scene>\n\n");			
			rsml.flush();

		}
		catch(IOException e){Util.log("CSV file creation failed ="+e);}
	}	
	
	
	/**
	 * Create CSV file which will contain the data from the simulation
	 * Create the headers.
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resequ : equivalent resistance
	 */
	public static void createCSV(boolean arch, boolean water, boolean aba, boolean carbon, boolean resequ, ParameterSet parms){
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
		Date d = new Date();
		date = dateFormat.format(d);
		
		try{
			planet = new PrintWriter(new FileWriter(parms.outputFileName));
			String header = "temps, hour, article, father, orgID, type, ordre, article_age, x, y, z, x1, y1, z1";
			if(arch) header = header.concat(",article_length, article_surface, article_volume, article_diameter, article_growth, article_dry_mass, dist_from_base");
			if(water) header = header.concat(",water_pot_endo, water_pot_exo, Kh, Lr, Lr_percent, radial_water_flux, axial_water_flux, aqp, cavitation");
			if(aba) header = header.concat(",aba_conc, aba_prod, aba_delta_conc, aba_delta_quant");
			if(carbon) header = header.concat(",maintenance_demand, growth_demand, dry_mass, growth_eff, photosynthesis");					  
			planet.println(header);
			planet.flush();
			if(resequ){
				resEqu = new PrintWriter(new FileWriter("/Users/guillaumelobet/Desktop/"+date+"_resequ.csv"));
				header = "temps, hour, article, type, resequ, water_pot_endo";
				resEqu.println(header);
				resEqu.flush();
			}   
			
		}
		catch(IOException e){Util.log("CSV file creation failed ="+e);}
	}
	

	  
	/**
	 * Update the csv file with the simulation data of the current time step
	 * @param node: list of plant nodes
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resequ: equivalent resistance
	 * @param t: additional time due to water resolution iterations (0-1)
	 */
	public static void RSMLExport(ArrayList<Node> node, float t){	
		rsml.println("	<plant>\n\n");	
		ArrayList<SegRoot> root1 = new ArrayList<SegRoot>(); 
		ArrayList<SegRoot> root2 = new ArrayList<SegRoot>(); 
		ArrayList<SegRoot> root3 = new ArrayList<SegRoot>(); 
		for (Node n : node) { 
			Article a = n.article;
			if(a instanceof SegRoot){
				SegRoot sr = (SegRoot) a;
				if(sr.getDistFromBase() == 0){
					if(sr.getOrdre() == 1) root1.add(sr);
					if(sr.getOrdre() == 2) root2.add(sr);
					if(sr.getOrdre() == 3) root3.add(sr);
				}
			}
		}
		for(SegRoot sr1 : root1){
			SegRoot a1 = sr1;
			while(a1.getParentArticle() != null){
				
				a1 = (SegRoot) a1.getParentArticle();
			}
			// Add the lateral roots
			for(SegRoot sr2 : root2){
				//if()
				SegRoot a2 = sr2;
				while(a2.getParentArticle() != null){
					
					a2 = (SegRoot) a2.getParentArticle();
				}
				
				root1.remove(a2);
			}
			
			root1.remove(a1);
		}
		
		Util.log("CSV export done in table "+date+"_planet.csv");
	}
	
	
	/**
	 * Update the csv file with the simulation data of the current time step
	 * @param node: list of plant nodes
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resequ: equivalent resistance
	 * @param t: additional time due to water resolution iterations (0-1)
	 */
	public static void CSVExport(ArrayList<Node> node, boolean arch, 
			boolean water, boolean aba, boolean carbon, boolean resequ, float t, ParameterSet parms){	  
		
		if(planet == null) createCSV(arch, water, aba, carbon, resequ, parms);
		
		Article p = null;
		
		for (Node n : node) { 
			Article a = n.article;
			if(a.getParentArticle() != null) p = a.getParentArticle();
			else p = a;
			
			String data = Time.getTime()+", "+(Time.getCurrentHour() + t)+", "+a.toString()+", "+p+", "+a.getArticleID()+", "+a.getType()+", "+a.getOrdre()+", "+a.getAge()+", "+
					a.node.x+", "+a.node.y+", "+a.node.z+", "+p.node.x+", "+p.node.y+", "+p.node.z;
			if(arch) data = data.concat(", "+(float) a.getLength()+", "+(float) a.getSurface()+", "+(float) a.getVolume()+", "+(float) a.getWidth()+", "+(float) a.getGrowth()+", "+
					  (float) a.getDryMass()+", "+(float) a.getDistFromBase());
			if(water) data = data.concat(", "+(float) a.getWaterPot()+", "+(float) a.getWaterPotExo()+", "+(float) a.getWaterAxialResistance()+", "+(float) a.getWaterRadialResistance()+", "+
					  (float) a.getLrPercent()+", "+(float) a.getRadialWaterFlux()+", "+(float) a.getAxialWaterFlux()+", "+ (float) a.getAQPCondEffect()+", "+ (float) a.getCavitationEffect());
			if(aba) data = data.concat(", "+(float) a.getConcSolute(1)+", "+(float) a.getABAProduction()+", "+(float) a.getDeltaConcSolute(1)+", "+(float) a.getDeltaQuantSolute(1));
			if(carbon) data = data.concat(", "+(float) a.getMaintenanceDemand()+", "+(float) a.getGrowthDemand()+", "+(float) a.getDryMass()+", "+(float) a.getGrowthEfficiency()+", "+(float) a.getPhotosynthesis());
			  
			planet.println(data);
			planet.flush();
		}
		//Util.log("CSV export done in "+parms.outputFileName);
		if(resequ){
			String dataRE = "";
			for (Node n : node) { 
				if (n.article.getParentArticle() instanceof ArtVirt){
					n.article.getEquivalentConductance();
					if(n.article instanceof SegRoot) dataRE = Time.getTime()+", "+Time.getCurrentHour()+", "+n.article+", 1, "+n.article.KhEqu+", "+n.article.getWaterPot();
					else dataRE = Time.getTime()+", "+Time.getCurrentHour()+", "+n.article+", 2, "+n.article.KhEqu+", "+n.article.getWaterPot();
					resEqu.println(dataRE);
					resEqu.flush();
				}
			}
		}
	}

	
	
	/**
	 * Create the SQL database
	 */
	public static void createDB(){
		try{
			Class.forName(driver).newInstance();
			java.sql.Connection conn = java.sql.DriverManager.getConnection("jdbc:mysql://localhost/mysql" ,usr, psw);
		  	java.sql.Statement stmt=conn.createStatement();
			stmt.executeUpdate("CREATE DATABASE "+database);
		}
		catch(Exception e){ System.out.println("Database is not created "+e);}	
		Util.log("Creation of the database "+database);
	}
	  
	/**
	 * Test if the SQL database exist. If not, it will be created
	 */
	public static void testDB(){
		  
		try{
			Class.forName(driver).newInstance();
			java.sql.DriverManager.getConnection(url ,usr, psw);
		}
		catch(Exception e){ createDB();}
	}
	  
	/**
	 * Create the SQL statement
	 * @return the statement
	 */
	public static Statement getSQLStatement(){
		  
		Connection con;
		Statement sql;
		  
		try { Class.forName(driver); }
		catch (ClassNotFoundException e) {
			Util.log("The driver " + driver + " was not found.");
			Util.log("You should check your Java installation.");
			Util.log("You will not be able to write to a database.");
			return null;
		}
		
		try {
			con = DriverManager.getConnection(url, usr, psw);
			sql = con.createStatement();
		}
		catch (SQLException sqlE) {
			Util.log("The specified database was not found.");
			Util.log("You will not be able to write to a database.");
			return null;
		}
		return sql;
	}

	/**
	 * Create the data table(s).
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resEq: equivalent resistance export
	 */
	public static void createDataTables(boolean arch, boolean water, boolean aba, boolean carbon, boolean resEq){
		  
		testDB();
		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmm");
		Date date = new Date();
		dbDate = (dateFormat.format(date));
		Statement sql = getSQLStatement();  	

		/*
		 * Table node_data 
		 */
		if(arch){
			tableND = "node_data_"+dbDate;
			try {
				sql.executeUpdate("CREATE TABLE " + tableND + " (temps INT);");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN hours INT;");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article CHAR(32);");
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN orgID INT;");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN father CHAR(32);");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN type INT;");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN ordre INT;");    
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_age FLOAT;");
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN x FLOAT;");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN y FLOAT;");  
				sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN z FLOAT;");
				if(arch){
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_length FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_volume FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_surface FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_diameter FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_growth FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN article_dry_mass FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN dist_from_base FLOAT;"); 
				}
				if(water){
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN water_pot_endo FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN water_pot_exo FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN Kh FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN Lr FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN Lr_percent FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN radial_water_flux DOUBLE;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN axial_water_flux DOUBLE;"); 
				}
				if(aba){
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN ABA_conc FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN ABA_prod FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN ABA_delta_quant FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN ABA_delta_conc FLOAT;");  
				}
				if(carbon){
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN maintenance_demand FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN growth_demand FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN dry_mass FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN growth_eff FLOAT;");  
					sql.executeUpdate("ALTER TABLE " + tableND + " ADD COLUMN photosynthesis FLOAT;");  
				}
			}
			catch (SQLException sqlE) {
				Util.log("The table " + tableND + " could not be created");
	            Util.log(sqlE.getMessage());
			}
			Util.log("Creation of the datatable "+tableND);
		}
		  
		/*
		 * Table res_equ 
		 */
		if(resEq){
			tableRE = "res_eq_"+dbDate;
			try {
				sql.executeUpdate("CREATE TABLE " + tableRE + " (temps INT, hour INT);");  
				sql.executeUpdate("ALTER TABLE " + tableRE + " ADD COLUMN article CHAR(32);");  
				sql.executeUpdate("ALTER TABLE " + tableRE + " ADD COLUMN type INT;");  
				sql.executeUpdate("ALTER TABLE " + tableRE + " ADD COLUMN resequ FLOAT;");  
				sql.executeUpdate("ALTER TABLE " + tableRE + " ADD COLUMN water_pot_endo FLOAT;");  
			}
			catch (SQLException sqlE) {
				Util.log("The table " + tableRE + " could not be created");
				Util.log(sqlE.getMessage());
			}
			Util.log("Creation of the datatable "+tableRE);
		}         
	  }

	
	/**
	 * Close the SQL connection
	 */
	public static void closeSQLConnection(){
		try{
			Statement sql = getSQLStatement();   
			sql.close();
		}
		catch(Exception e){	Util.log("Connection closure failed in Exports.java"); }
	}
	  
	
	/**
	 * Export the parameters of the simulation
	 */
	public static void exportSimulationParameters(ParameterSet p){
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		String day = dateFormat.format(new Date());
		  
		DateFormat hourFormat = new SimpleDateFormat("HH:mm");
		String hour = hourFormat.format(new Date());
		  
		DateFormat datehourFormat = new SimpleDateFormat("yyyyMMdd_HHmm");
		dbDate = (datehourFormat.format(new Date()));
		tableND = "node_data_"+dbDate;
		
		String request1 = ""; 
		  
		try{
			Statement sql = getSQLStatement();
			request1 = "INSERT INTO "+tableParam+" ";
			String request2 = " (table_name, date, hour, base_Temp, time_step, export_type, export_noon, export_arch, export_water, export_aba, export_carbon, "+
					  			" export_res_eq, start_Export, resolve_water, hydraulic_lift, resolve_water_type, resolve_iteration, rewatering_frequence, aquaporin, "+
					  			" cavitation, stomata, stress_type, start_stress, radial_mod, axial_mod, aba_solve_method, aba_production, resolve_carbon, variable_partitionning, "+
					  			" root_shoot_partitionning, root_tortuosity, diameter_stochasticity, branching_level, seminal_number, adventive_roots, seminal_angle_stochasticity, "+
					  			" primary_max_length, secondary_max_length, tertiary_max_length, primary_lauz, secondary_lauz, primary_interbranch, secondary_interbranch, "+
					  			" secondary_insert_angle, tertiary_insert_angle, primary_init_diameter, seminal_init_diameter, crown_init_diameter, principal_gravitropism, "+
					  			" seminal_gravitropism, secondary_gravitropism, tertiary_gravitropism)";
			String request3 = " VALUES ('"+tableND+"', '"+day+"', '"+hour+"', "+p.baseTemp+", "+p.timeStep+", "+p.exportType+", "+p.exportAtNoon+", "+p.exportArchitecture+
					  			", "+p.exportWater+", "+p.exportABA+", "+p.exportCarbon+", "+p.exportResEq+", "+p.startExport+", "+p.resolveWater+", "+p.allowHydraulicLift+
					  			", "+p.resolveWaterType+", "+p.resolveIteration+", "+p.reWateringFrequence+", "+p.aqp+", "+p.cavitation+", "+p.stomata+", "+p.stressType+", "+p.startStress+
					  			", "+p.radMod+", "+p.axMod+", "+p.ABASolveMethod+", "+p.ABAProduction+", "+p.resolveCarbon+", "+p.variablePartitionning+", "+p.rootShootPartitionning+", "+p.rootTortuosity+
					  			", "+p.diamStochasticity+", "+p.branchingLevel+", "+p.numberOfSeminals+", "+p.adventives+", "+p.semAngStochasticity+", "+p.primaryMaxLength+", "+p.secondaryMaxLength+
					  			", "+p.tertiaryMaxLength+", "+p.primaryLAUZ+", "+p.secondaryLAUZ+", "+p.primaryInterBranch+", "+p.secondaryInterBranch+", "+p.secondaryInsertAngle+", "+p.tertiaryInsertAngle+
					  			", "+p.principalRootInitDiam+", "+p.seminalRootInitDiam+", "+p.crownRootInitDiam+", "+p.principalRootGravitropism+", "+p.seminalRootGravitropism+", "+p.secondaryRootGravitropism+
					  			", "+p.tertiaryRootGravitropism+")";
		  		
			request1 = request1.concat(request2.concat(request3));
			sql.execute(request1);
			  
		}
		catch(Exception e){
			Util.log("Export simulation parameters faileld");
			Util.log(""+e);
			Util.log(""+request1);
		}				
	}
	
	
	/**
	 * Update the "simulation_duration" parameters in the simulation table
	 */
	public static void updateExportSimulationParameters(){
		String request1 = "";
		try{
			Statement sql = getSQLStatement();
			request1 = "UPDATE  "+tableParam+" SET";
			String request2 = " simulation_duration = "+Time.getTime();
			String request3 = " WHERE table_name = '"+tableND+"'";
			  		
			request1 = request1.concat(request2.concat(request3));
			sql.execute(request1);
			  
		}
		catch(Exception e){
			Util.log("Finalize export simulation parameters faileld");
			Util.log(""+e);
			Util.log(""+request1);
		}		
	}
	
	/**
	 * Export the simulation data to the SQL table.
	 * @param node: list of plant nodes
	 * @param arch: architectural data
	 * @param water: water related data
	 * @param aba: aba related data
	 * @param carbon: carbon related data
	 * @param resEq: equivalent resistance export
	 */
	  public static void SQLexport(ArrayList<Node> node, boolean arch, boolean water, boolean aba, boolean carbon, boolean resEq){
		  
		  DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd hh:mm");
		  Date d = new Date();
		  date = dateFormat.format(d);
		  
		  ArrayList<Node> nodeList = node;			  		
		  Article a = null;
		  MerisRoot mr = null;
		  double radialRes = 0;
		  double axRes = 0;
		  String f = "";

		  /*
		   * Table node_data
		   */		  
		  if(arch){
			  try{
				  Statement sql = getSQLStatement();
				  if(arch){
					  for (Node n : nodeList) { 
						  if(n.article instanceof MerisRoot){
							  mr = (MerisRoot) n.article;
							  mr.setDistanceFromBase();
						  }
					  }
				  }
				  for (Node n : nodeList) { 
					  a = n.article;				  		
					  if(a.getParentArticle() != null) f = a.getParentArticle().toString();
					  else f = "";
				  		
					  String request1 = "INSERT INTO "+tableND+" ";
					  String request2 = " (article, father, orgID, type, ordre, temps, hours, article_age, x, y, z";
					  String request3 = " VALUES ('"+a.toString()+"', '"+f+"', "+a.getArticleID()+", "+a.getType()+", "+a.getOrdre()+
							  ", "+Time.getTime()+", "+Time.getCurrentHour()+", "+a.getAge()+", "+a.node.x+", "+a.node.y+", "+a.node.z;
				  		
					  if(arch) {
						  request2 = request2.concat(",article_length, article_surface, article_volume, article_diameter, article_growth, article_dry_mass, dist_from_base");
						  request3 = request3.concat(", "+a.getLength()+", "+a.getSurface()+", "+a.getVolume()+", "+a.getWidth()+", "+a.getGrowth()+
								  ", "+a.getDryMass()+", "+a.getDistFromBase());
					  }
					  if(water){
						  if(1/a.getWaterRadialResistance() == 0) radialRes = 1e10;
						  else radialRes = a.getWaterRadialResistance();				  		
						  if(1/a.getWaterAxialResistance() == 0) axRes = 1e10;
						  else axRes = a.getWaterAxialResistance();
						  request2 = request2.concat(",water_pot_endo, water_pot_exo, Kh, Lr, Lr_percent, radial_water_flux, axial_water_flux");
						  request3 = request3.concat(", "+a.getWaterPot()+", "+a.getWaterPotExo()+", "+axRes+", "+radialRes+", "+a.getLrPercent()+", "+a.getRadialWaterFlux()+", "+a.getAxialWaterFlux());
					  }
					  if(aba){
						  request2 = request2.concat(",aba_conc, aba_prod, aba_delta_conc, aba_delta_quant");
						  request3 = request3.concat(", "+a.getABAConsumption()+", "+a.getABAProduction()+", "+a.getDeltaConcSolute(1)+", "+a.getDeltaQuantSolute(1));
					  }
					  if(carbon){
						  request2 = request2.concat(",maintenance_demand, growth_demand, dry_mass, growth_eff, photosynthesis");
						  request3 = request3.concat(", "+(float) a.getMaintenanceDemand()+", "+(float) a.getGrowthDemand()+", "+(float) a.getDryMass()+", "+(float) a.getGrowthEfficiency()+", "+(float) a.getPhotosynthesis());
					  }
					  request2 = request2.concat(")");
					  request3 = request3.concat(")");
				  		
					  request1 = request1.concat(request2.concat(request3));
					  sql.execute(request1);
				  }
			  }
			  catch(Exception e){ Util.log("Export to "+tableND+" failed in Exports.java: "+e);	}
			  Util.log("Data export for time "+Time.getTime()+" done in table "+tableND);
		  }
		  
		  /*
		   * Table resequ
		   */
		  if(resEq){
			  try{
				  Statement sql = getSQLStatement();
				  			    
				  for (Node n : nodeList) { 
					  if (n.article.getParentArticle() instanceof ArtVirt){
						  n.article.getEquivalentConductance();
						  if(n.article instanceof SegRoot){
							  sql.execute("INSERT INTO "+tableRE+" (temps, hour, article, type, resequ, water_pot_endo) " +
									  "VALUES ("+Time.getTime()+", "+Time.getCurrentHour()+", '"+n.article+"', 1, "+n.article.KhEqu+", "+n.article.getWaterPot()+" )");
						  }
						  else{
							  sql.execute("INSERT INTO "+tableRE+" (temps, hour, article, type, resequ, water_pot_endo) " +
									  "VALUES ("+Time.getTime()+", "+Time.getCurrentHour()+", '"+n.article+"', 2, "+n.article.KhEqu+", "+n.article.getWaterPot()+" )");
						  }
					  }
				  }
			  }
			  catch(Exception e){ Util.log("Export to "+tableRE+" failed in Exports.java: "+e); }
			  Util.log("Data export for time "+Time.getTime()+" done in table "+tableRE);
		  }
		  closeSQLConnection();
	  }
	}
