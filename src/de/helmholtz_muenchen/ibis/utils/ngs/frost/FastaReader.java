package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * @author tanzeem.haque
 *
 */


public final class FastaReader{

    /**
     * After much Ueberlegung each sequence is doch represented in a String
     */
	
	
    private String identifier;
    private String sequence;
    
//    private String fileName;
//    private ArrayList<LinkedList<Character>> sequence;

    
    public FastaReader() {
    	
		// TODO Auto-generated constructor stub

	}
    
    void readSequenceFromFile(String file, String currentChr) {
		
//		File file = new File(fastaFile);

		String id = "";
		String seq = "";


		
		try { 
			
			/**File input = new File(file);
			Scanner sc = new Scanner(input, "UTF-8");**/

			BufferedReader in = new BufferedReader(new FileReader(file), 10000000);
				        String currentLine = "";

	        StringBuffer buffer = new StringBuffer();
	        /*
	        if((currentLine = in.readLine()) == null )
	            throw new IOException( file + " File is empty " );
	     
	        if(currentLine.charAt(0) != '>')
	            throw new IOException("Header of the file " + file + " must be '>'");
	        else {
	        	if (currentChr.equals(currentLine.substring(1).trim()))//.replaceAll("\\s+", "_"))) 
	        		{
	        		id=currentChr;
	        		}
//	        	System.err.println(currentLine.substring(1).trim());}
	        }
	        */

        	boolean read = false;
	        for(currentLine = in.readLine().trim(); currentLine != null; currentLine = in.readLine()) {
	        	/** while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();**/
//	        	String tmp_id = "";//currentLine.substring(1).trim();//.replaceAll("\\s+", "_");
//	        	System.err.println(tmp_id);
	        	if(currentLine.length() > 0 && currentLine.charAt(0) == '>') {
	        		if (currentLine.substring(1).trim()./*replaceAll("\\s+", "_").*/equals(currentChr)) {
						seq = buffer.toString();
	        			buffer = new StringBuffer();
		        		id=currentLine.substring(1).trim();
//		        		System.out.println("ID in FR: " + id);
		        		read = true;
//		        		FrostRunner.memory();
	        		}
	        		else {
	        			read = false;
//	        			System.out.println("NOT READING");
		        		continue;
	        		}
	        	}
	        		
	        	else {
					if (read) {
						buffer.append(currentLine.trim());

					}
//					else
//						continue;
	        		
        		}
	        }
//    		FrostRunner.memory();

		
	        if( buffer.length() != 0 ){
        		seq = buffer.toString(); 
//        		System.out.println("Sequence: " + seq.substring(seq.length()-10, seq.length()-1));
	        }

		}
		catch(IOException e){ 
			System.out.println("Error when reading "+file);
	        e.printStackTrace();
	        
		}
		identifier = (String) id;
		sequence = (String) seq;
		
	}
	
//    void readSequenceFromFile(String file) {
    	//
//    			String id = "";
//    			String seq = "";
//    	    	try { 
//    	    		BufferedReader in = new BufferedReader(new FileReader(file));
//    	    	    StringBuffer buffer = new StringBuffer();
//    	    	    String currentLine = "";
//    	    	    if((currentLine = in.readLine()) == null )
//    	    	        throw new IOException( file + " File is empty " );
//    	    	 
////    	    	    if(currentLine.charAt(0) != '>')
////    	    	        throw new IOException("Header of the file " + file + " must be '>'");
////    	    	    else
//    	    	        id = currentLine.substring(1).trim(); //.replaceAll("\\s+", "_");
//    	    	    
    	//
//    	    	    for(currentLine = in.readLine().trim(); currentLine != null; currentLine = in.readLine()) {
//    	    	    	if(currentLine.length() > 0 && currentLine.charAt(0) == '>'){
//    	    	    		seq = buffer.toString();
//    	    	    		buffer = new StringBuffer();
//    	    	    		id = currentLine.substring(1).trim(); //.replaceAll("\\s+", "_");
//    	    	    		
//    	    	    	} 
//    	    	    	else {
//    	    	    		buffer.append(currentLine.trim());
//    	    	    	}
//    	    	    	
//    	    	    }
    	//
//    	    	    if( buffer.length() != 0 ){
//    	    			seq = buffer.toString();     		
//    	    	    }
    	//
//    	    	}
//    	    	catch(IOException e){ 
//    	    		System.out.println("Error when reading "+file);
//    	    	    e.printStackTrace();
//    	    	    
//    	    	}
//    	    	identifier = (String) id;
//    			sequence = (String) seq;
//    	    	
//    	    }
    /**
	 * @return first sequence as a char linked list
	 */
	
	String getSequence() { 
		return sequence;
		
	}

	/**
	 * @return first identifier as String

	 */
	String getIdentifier() {
		return identifier;
		
	}	

	/**
	 * @return the length of the first sequence
	 */
	int getLength() {
		return sequence.length();
	}
	

	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return ">" + getIdentifier() + "\n" + getSequence() + "\n";
	}

	
	
}



/**
 * This class will read first sequence from a Fasta format file 
 */









//public final class FastaReader{
//
//    /**
//     * After much Ueberlegung each sequence is doch represented in a String
//     */
//	
//	
//    private String identifier;
//    private String sequence;
//    
////    private String fileName;
////    private ArrayList<LinkedList<Character>> sequence;
//
//    
//    public FastaReader() {
//    	
//		// TODO Auto-generated constructor stub
//
//	}
//    
//    void readSequenceFromFile(String fileName, String currentChr) {
//
//		String id = "";
//		String seq = "";
//
//		try { 
//    		BufferedReader in = new BufferedReader(new FileReader(fileName));
//    	    StringBuffer buffer = new StringBuffer(50);
//    	    String currentLine = "";
//		    Boolean read = false;
//
////    	    if((currentLine = in.readLine()) == null )
////    	        throw new IOException( fileName + " File is empty " );
////    	 
////    	    if(currentLine.charAt(0) != '>')
////    	        throw new IOException("Header of the file " + fileName + " must be '>'");
//
//    	    while((currentLine = in.readLine()) != null) {
//    	    	if(currentLine.length() > 0 && currentLine.charAt(0) == '>'){
////    	    		System.out.println("I'M READING HEADER");
//    	    		if(currentLine.substring(1).equals(currentChr)){
//		    			seq = buffer.toString();
//		    			buffer = new StringBuffer();
//		    			id = currentLine.substring(1).trim(); //.replaceAll("\\s+", "_");
////		    			System.out.println("ID FASTA: " + id);
//		    			read = true;
//		    		}
//		    		else 
//		    			read = false; 
//    	    	}
//    	    	if(read) {
//    	    		buffer.append(currentLine.trim());
//    	    	} 
//    	    	
//    	    	else 
//    	    		break;
//    	    }
//    	    if( buffer.length() != 0 ){
//    			seq = buffer.toString();     		
//    	    }
//		}
//		catch(IOException e){ 
//			System.out.println("Error when reading "+fileName);
//	        e.printStackTrace();
//	        
//		}
//		identifier = (String) id;
//		sequence = (String) seq;
//		
//	}
//	
//    void readSequenceFromFile(String file) {
//
//		String id = "";
//		String seq = "";
//    	try { 
//    		BufferedReader in = new BufferedReader(new FileReader(file));
//    	    StringBuffer buffer = new StringBuffer();
//    	    String currentLine = "";
//    	    if((currentLine = in.readLine()) == null )
//    	        throw new IOException( file + " File is empty " );
//    	 
////    	    if(currentLine.charAt(0) != '>')
////    	        throw new IOException("Header of the file " + file + " must be '>'");
////    	    else
//    	        id = currentLine.substring(1).trim(); //.replaceAll("\\s+", "_");
//    	    
//
//    	    for(currentLine = in.readLine().trim(); currentLine != null; currentLine = in.readLine()) {
//    	    	if(currentLine.length() > 0 && currentLine.charAt(0) == '>'){
//    	    		seq = buffer.toString();
//    	    		buffer = new StringBuffer();
//    	    		id = currentLine.substring(1).trim(); //.replaceAll("\\s+", "_");
//    	    		
//    	    	} 
//    	    	else {
//    	    		buffer.append(currentLine.trim());
//    	    	}
//    	    	
//    	    }
//
//    	    if( buffer.length() != 0 ){
//    			seq = buffer.toString();     		
//    	    }
//
//    	}
//    	catch(IOException e){ 
//    		System.out.println("Error when reading "+file);
//    	    e.printStackTrace();
//    	    
//    	}
//    	identifier = (String) id;
//		sequence = (String) seq;
//    	
//    }
//	/**
//	 * @return first sequence as a char linked list
//	 */
//	
//	String getSequence() { 
//		return sequence;
//		
//	}
//
//	/**
//	 * @return first identifier as String
//
//	 */
//	String getIdentifier() {
//		return identifier;
//		
//	}	
//
//	/**
//	 * @return the length of the first sequence
//	 */
//	int getLength() {
//		return sequence.length();
//	}
//	
//
//	
//	/* (non-Javadoc)
//	 * @see java.lang.Object#toString()
//	 */
//	@Override
//	public String toString() {
//		// TODO Auto-generated method stub
//		return ">" + getIdentifier() + "\n" + getSequence() + "\n";
//	}
//
//	
//}



/**
 * This class will read first sequence from a Fasta format file 
 */
