/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
