package de.helmholtz_muenchen.ibis.utils.ngs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Pattern;


public class FileValidator {
	
	public static final String SAMTOOLS_VERSION = "0.1.18";
	public static final String BWA_VERSION = "0.7.5a";
	public static final String BFAST_VERSION = "0.7.0";
	
	
	
	public FileValidator(String path2file){

	}
	
	
	/**
	 * Verifies fasta file format
	 * @return true if file is in fasta format and is for nucleotide sequences
	 */
	public static boolean checkFastaFormat(String path2file){
		
		BufferedReader br = null;
		try {

			String sCurrentLine;
			br = new BufferedReader(new FileReader(path2file));
			boolean stop = false;
			
			while (!stop && (sCurrentLine = br.readLine())!= null) {
				if(!(sCurrentLine.equals(""))){//if not an empty line
					boolean header = (sCurrentLine.charAt(0)=='>');
					sCurrentLine = br.readLine();
					boolean seq = (Pattern.matches("^([A,C,G,T,U,R,Y,K,M,S,W,B,D,H,V,N,X])+$", sCurrentLine.toUpperCase()));
					if(header & seq){
						//System.out.println("Inputfile IS in fasta format");
						return true;
					}
					stop=true; 	
				}
			}		
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return false;
	}
	
	/**
	 * Checks, if a file is in SAM or BAM format
	 * @param path2File
	 * @return
	 */
	public static boolean checkSamBamFormat(String path2File){
		// TODO: implement check
		return true;
	}
	
	
	/**
	 * Verifies fastq file format by checking only the first (!) entry
	 * @return true if file is in fastq format and is for nucleotide sequences
	 */
	public static boolean checkFastqFormat(String fastqPath){
		BufferedReader br = null;
		try {

			String sCurrentLine;
			br = new BufferedReader(new FileReader(fastqPath));
			boolean stop = false;
			
			while (!stop && (sCurrentLine = br.readLine())!= null) {
				if(!(sCurrentLine.equals(""))){//if not an empty line
					//boolean header = (Pattern.matches("^@[A-Za-z0-9_\\.:\\-\\/\\s]*", sCurrentLine));
					boolean header =  (sCurrentLine.charAt(0)=='@');
					sCurrentLine = br.readLine();
					boolean seq = (Pattern.matches("^[A,C,G,T,U,R,Y,K,M,S,W,B,D,H,V,N,X]+$", sCurrentLine.toUpperCase()));
					sCurrentLine = br.readLine();
					//boolean seqname = (Pattern.matches("^\\+[A-Za-z0-9_\\.:\\-\\/\\s]*", sCurrentLine));
					boolean seqname =  (sCurrentLine.charAt(0)=='+');
					sCurrentLine = br.readLine();
					boolean qual = (Pattern.matches("^[\\S]+$", sCurrentLine));
					if(header && seqname && seq && qual){
						return true;
					}else{
						System.out.println("Header is correct: "+header);
						System.out.println("Seqname is correct: "+seqname);
						System.out.println("Sequence is correct: "+seq);
						System.out.println("Quality values are correct: "+qual);
					}
					stop=true; 
				}
			}		
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return false;
	}
	
	/**
	 * Verifies the barcode file needed as input for RawReadManipulator
	 * @return true if file has right format TODO: At the moment only a rough validation is performed. Maybe improve it! Check if header line?
	 */
	public static boolean validateBarcodeFile(String filepath){
		
		BufferedReader br = null;
		try {

			String sCurrentLine;
			br = new BufferedReader(new FileReader(filepath));
			
			while ((sCurrentLine = br.readLine())!= null) {
				if(!(sCurrentLine.equals(""))){//if not an empty line
					boolean does_match = (Pattern.matches("^\\w+\\s\\w+$", sCurrentLine));
					if(!does_match){
						return false;
					}

				}
			}
			return true;
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return false;
	}

	/**
	 * Checks the available user version of the input program and compares it to our version
	 * @param path2program and name of the tool(SAMTOOLS,BWA,BFAST)
	 * @return 0 if Version are equal, 1 if user has a new version, 2 if user has older version
	 */
	public static int versionControl(String path2tools,String toolname){
        String com = path2tools;
        String VERSION="";
        if(toolname.equals("BWA")){
            VERSION=BWA_VERSION;	
        }else if(toolname.equals("SAMTOOLS")){
        	 VERSION=SAMTOOLS_VERSION;	
        }else if(toolname.equals("BFAST")){
        	 VERSION=BFAST_VERSION;	
        }
        
        Process p;
		try {
			p = Runtime.getRuntime().exec(com);
	        p.waitFor();
	      	BufferedReader stdError = new BufferedReader(new 
	                  InputStreamReader(p.getErrorStream()));
	      	String s = null;
	          
	          
	         while ((s = stdError.readLine()) != null) {
	        	 //Find the version line
	              boolean line = Pattern.matches("(?m)^Version.*$", s);
	              if(line){
	    	          String[] fields =  s.split("\\s");
	    	          String curr_version = fields[1];
	    	          
	    	          //if versions are equal return 0
	    	          if(VERSION.equals(curr_version)){
	    	        	 return 0;
	    	          }else{// else compare
	    	        	  String[] v1=VERSION.split("\\.");
	    	        	  String[] v2={""};  
	    	              if(toolname.equals("BWA")){
	    	            	  String[]f = curr_version.split("-");
	    	            	  if(f[0].equals(VERSION)){
	    	            		  return 0;
	    	            	  }
	    	            	  v1 = VERSION.split("\\.");
	    	            	  v2 = f[0].split("\\.");	
	    	              }else if(toolname.equals("SAMTOOLS")){
	    	            	  String[] v_temp = curr_version.split("-");
	    	            	  v2 = v_temp[0].split("\\.");
	    	              }else if(toolname.equals("BFAST")){
	    	            	  char last = curr_version.charAt(curr_version.length()-1);
	    	            	  boolean isAnInt= last>='0' && last<='9';
	    	            	  if(!isAnInt){
	    	            		  curr_version=curr_version.substring(0, curr_version.length()-1);  
	    	            	  }
	    	            	  if(VERSION.equals(curr_version)){
	    	            		  return 0;
	    	            	  }
	    	            	  v2 = curr_version.split("\\.");
	    	              }else{
	    	            	  v1=null;
	    	            	  v2=null;
	    	              }


	    	              
	    	        	  for(int i=0;i< v1.length;i++){
	    	        		  
		    	              //Remove all non-digits from the strings
	    	        		  // This means the version control is not exact!
	    	        		  Integer n1 = Integer.parseInt(v1[i].replaceAll("\\D*", ""));
	    	        		  Integer n2 = Integer.parseInt(v2[i].replaceAll("\\D*", ""));
		    	              
	    	        		  if(n1 > n2){ //User version is older 
	    	        			  return 2;
	    	        		  }else if(n1 < n2){//User version is newer
	    	        			  return 1;
	    	        		  }
	    	        	  }
	    	        	 
	    	          }
	              }
	          }       
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}


		return -1;
	}
	
	
}
