package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import de.helmholtz_muenchen.ibis.ngs.FrostslilheplerNodeModel;


public class Chunker {

	/**
	 * @param args
	 * This class will split ONE single fasta according to a given chunk length into multiple chunks.
	 * Hekper class for FrostsL'lHelper
	 */
	private String file;
	
	public Chunker(String file) {
		setFile(file);
	}	
	public String getFile() {
		return file;
	}
	public void setFile(String file) {
		this.file = file;
	}
	
	public void createChunk() throws IOException {
		String s[] = this.file.split("/");
		
		FastaReader fr = new FastaReader();
		String id = getID(this.file);
//		System.out.println(id);
		fr.readSequenceFromFile(this.file, id);
		int seqLength = fr.getLength();
//		System.out.println(seqLength);
		int chunk = (seqLength/FrostslilheplerNodeModel.chunk_length) +1;
//		System.out.println(chunk);
		
		for (int i = 0; i < chunk; i++) {
			
			int begin = i*FrostslilheplerNodeModel.chunk_length;
			int end = begin + FrostslilheplerNodeModel.chunk_length;
			if (end > seqLength)
				end = seqLength;
			
//			System.out.println("Chunk :" + i + "\t" + "Begin: " + begin + "\t" + "End: " + end);
			String tmp = i + "_" + s[s.length-1];
			String newFile = this.file.replace(s[s.length-1], tmp);
			
			BufferedWriter bw= new BufferedWriter(new FileWriter(newFile), 10000000);
			try { 
//				System.out.println(newFile);
				bw.write(fr.getSequence().substring(begin, end));
			} catch (Exception e) {
				System.err.println("Error: " + e.getMessage());
				e.printStackTrace();
			}
			bw.close();
		}
		
	}
	
	private String getID(String file) {
		
		String id = "";
		
		try { 

			BufferedReader in = new BufferedReader(new FileReader(file), 10000000);
				        String currentLine = "";
	       
	        for(currentLine = in.readLine().trim(); currentLine != null; currentLine = in.readLine()) {

	        	if(currentLine.length() > 0 && currentLine.charAt(0) == '>') {
		        	id=currentLine.substring(1).trim();
	        		
	        	}	        		
	        	
	        }
		}
		catch(IOException e){ 
			System.out.println("Error when reading "+file);
	        e.printStackTrace();
	        
		}
		return id;
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		Chunker chunker = new Chunker("/home/ibis/tanzeem.haque/Documents/Frost_outputs/chr21_C0.fa");
		chunker.createChunk();

	}

}
