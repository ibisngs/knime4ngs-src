package de.helmholtz_muenchen.ibis.utils.ngs;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;


public class PEDFile {

	
	private HashMap<String, String> PED_FILE = new HashMap<String, String>();
	
	public PEDFile(String PATH_TO_FILE){
		
		try{
  		  FileInputStream fstream = new FileInputStream(PATH_TO_FILE);
  		  // Get the object of DataInputStream
  		  DataInputStream in = new DataInputStream(fstream);
  		  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  		  String strLine;
  		  //Read File Line By Line
  		  
//			     Family ID
//			     Individual ID
//			     Paternal ID
//			     Maternal ID
//			     Sex (1=male; 2=female; other=unknown)
//			     Phenotype

  		  while ((strLine = br.readLine()) != null)   {
  			  String[] fields = strLine.split("\t");
  			  
  			  String IndID = fields[1];
  			  String PaternalID = fields[2];
  			  String MaternalID = fields[3];
  			  String Phenotype = fields[5];
  			  
  			  if(Phenotype.equals("2")){				//2 == Phenotype 'affected'
  				  PED_FILE.put(IndID, PaternalID+"\t"+MaternalID);
  			  }		  
  		  }
  		  //Close the input stream
  		  in.close();
  		}catch (Exception e){//Catch exception if any
  			System.err.println("Error: " + e.getMessage());
  		}
	}
	
	public HashMap<String,String> getPED(){
		return PED_FILE;
	}
	
}
