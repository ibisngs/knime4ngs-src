package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * @author tanzeem.haque
 *
 */
public final class FastaReader{

    /**
     * After much Ueberlegung each sequence is doch represented in a String
     */
	
	
    private String [] identifier;
    private String [] sequence;
    
//    private String fileName;
//    private ArrayList<LinkedList<Character>> sequence;

    
    public FastaReader() {
    	
		// TODO Auto-generated constructor stub

	}
    
    void readSequenceFromFile(String file) {
		
//		File file = new File(fastaFile);

		ArrayList<String> id = new ArrayList<String>();
		ArrayList<String> seq = new ArrayList<String>();


		
		try { 
			BufferedReader in = new BufferedReader(new FileReader(file));
	        StringBuffer buffer = new StringBuffer();
	        String currentLine = "";
	        if((currentLine = in.readLine()) == null )
	            throw new IOException( file + " File is empty " );
	     
	        if(currentLine.charAt(0) != '>')
	            throw new IOException("Header of the file " + file + " must be '>'");
	        else
	            id.add(currentLine.substring(1).trim().replaceAll("\\s+", "_"));
	        

	        for(currentLine = in.readLine().trim(); currentLine != null; currentLine = in.readLine()) {
	        	if(currentLine.length() > 0 && currentLine.charAt(0) == '>'){
	        		String tmp_s = buffer.toString();
	        		seq.add(tmp_s);
	        		buffer = new StringBuffer();
	        		id.add(currentLine.substring(1).trim().replaceAll("\\s+", "_"));
	        		
	        	} 
	        	else {
	        		buffer.append(currentLine.trim());
	        	}
	        	
	        }
		
	        if( buffer.length() != 0 ){
        		String tmp_s = buffer.toString();     		
        		seq.add(tmp_s);

	        }

		}
		catch(IOException e){ 
			System.out.println("Error when reading "+file);
	        e.printStackTrace();
	        
		}
		//seq.size() = #sequences
		identifier = new String[id.size()];
		sequence = new String[seq.size()];

		for (int i = 0; i < seq.size(); i++) {
			identifier[i] = (String) id.get(i);
			sequence[i] = (String) seq.get(i);
		}
	}
	
	/**
	 * @return first sequence as a char linked list
	 */
	
	String getSequence() { 
		return sequence[0];
		
	}

	/**
	 * @return first identifier as String

	 */
	String getIdentifier() {
		return identifier[0];
		
	}	

	/**
	 * @return the length of the first sequence
	 */
	int getLength() {
		return sequence[0].length();
	}
	
	/**
	 * @param i
	 * @return sequence as a char linked list

	 */
	String getSequence(int i) { 
		return sequence[i];
		
	}
	
	/**
	 * @param i
	 * @return identifier as String

	 */
	String getIdentifier(int i) {
		return identifier[i];
		
	}
	
	/**
	 * @param i
	 * @return the length of the sequence at index i
	 */
	int getLength(int i) {
		return sequence[i].length();
	}
	/**
	 * @return number of sequences in the multifasta
	 */
	int size() {
		return sequence.length;
		
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		String data = "";
		System.out.println(size());
		System.out.println(identifier.length);
		for (int i = 0; i < size(); i++) {
			data += getIdentifier(i) + "\n";
			for (int j = 0; j < getSequence(i).length(); j++) {
				data += getSequence(i) + " ";
			}
			data += "\n";
		}
		return data;
	}

	
}



/**
 * This class will read first sequence from a Fasta format file 
 */
