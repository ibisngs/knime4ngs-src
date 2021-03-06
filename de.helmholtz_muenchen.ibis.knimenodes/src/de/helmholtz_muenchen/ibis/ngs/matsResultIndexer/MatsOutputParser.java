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


package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Abstract class which reads the MATS annotation and can convert it to gff3 which MISO can read.
 * @author michael.kluge
 *
 */
public abstract class MatsOutputParser {
	
	protected static final String NAME_FDR = "FDR";
	public static final String GENNAME = ";Genename=";
	public static final String COMMENT_SEP = ";";
	public static final String NAME = "Name=";
	public static final String NEWLINE = "\n";
	public static final String TAB = "\t";
	private static final String SEP_NAME = "@";
	private static final String GENE = "gene";
	protected final static HashMap<String, ArrayList<String>> NAMES = new HashMap<String, ArrayList<String>>(); 
	protected final HashMap<String, Integer> NAME2ID = new HashMap<String, Integer>();
	protected final StringBuffer OUT = new StringBuffer();
	
	private BufferedReader BR;
	protected BufferedWriter BW;
	protected boolean SILENT_MODE = false;
	protected final Double FDR;
	protected final ArrayList<String> EVENT_NAMES = new ArrayList<String>();
	
	/**
	 * Constructor
	 * @param path
	 * @throws IOException
	 */
	public MatsOutputParser(String path, BufferedWriter bw, Double FDR) throws IOException {
		File f = this.openFile(path);
		
		if(f != null) {
			this.BR = new BufferedReader(new FileReader(f));
			this.BW = bw;
			
			// check, if silence mode is on
			if(FDR != null && FDR <= 1 & FDR >= 0) {
				this.SILENT_MODE = true;
				this.FDR = FDR;
			}
			else 
				this.FDR = null;
			
			this.processFile();
		}
		else 
			this.FDR = null;
	}
	
	/**
	 * opens the file, if possible
	 * @param path
	 * @return
	 */
	protected File openFile(String path) {
		File f = new File(path);
		if(f.exists() & f.canRead()) {
			return f;
		}
		return null;
	}
	
	/**
	 * processes the file line by line
	 * @throws IOException
	 */
	protected void processFile() throws IOException {
		if(this.BR != null) {
			String head = this.BR.readLine();
			
			int length = this.init(head);
			
			String line;
			while((line = this.BR.readLine()) != null) {
				String split[] = line.split(TAB);
				
				if(split.length == length) {
					String name = this.getEvent(split);
					
					// add the name, if in silent mode
					if(this.SILENT_MODE && name != null)
						EVENT_NAMES.add(name);
				}
				else 
					throw new IllegalArgumentException("Line '" + line + "' has not enough or too much cols.");
			}
			this.closeFile();
		}
	}
	
	/**
	 * is used to find the right col names if annotation might ever change
	 * @param head
	 * @return
	 */
	protected int init(String head) {
		String split[] = head.split(TAB);
		boolean found;
		// get names to find
		ArrayList<String> names = NAMES.get(this.getType());
		if(this.SILENT_MODE)
			names.add(NAME_FDR);
		
		for(String n : names) {
			found = false;
			for(int i = 0; i < split.length; i++) {
				if(split[i].equals(n)) {
					this.NAME2ID.put(n, i);
					found = true;
				}
			}
			
			if(!found) {
				throw new IllegalArgumentException("Col '" + n + "' was not found in header line of file.");
			}
		}
		return split.length;
	}
	
	/**
	 * closes the file handle
	 */
	protected void closeFile() {
		if(this.BR != null) {
			try {
				this.BR.close();
			} catch (IOException e) { e.printStackTrace(); }
			this.BR = null;
		}
	}
	
	/**
	 * GFF format of the file
	 * @return
	 */
	public String getGFF() {
		return this.OUT.toString();
	}
	
	/**
	 * creates the gene line of the gff file for a gene
	 * @param gene
	 * @param chr
	 * @param exon
	 * @return
	 */
	protected String getGeneLine(String gene, String chr, Exon ...exon) {
		StringBuffer b = new StringBuffer();
		b.append(chr);
		b.append(TAB);
		b.append(this.getType());
		b.append(TAB);
		b.append(GENE);
		b.append(TAB);
		b.append(exon[0].getStart());
		b.append(TAB);
		b.append((exon[exon.length - 1]).getEnd());
		b.append(TAB);
		b.append(Isoform.SPACER);
		b.append(TAB);
		b.append(exon[0].getStrand());
		b.append(TAB);
		b.append(Isoform.SPACER);
		b.append(TAB);
		
		// add comment stuff;
		String name = this.getName(exon);
		b.append(NAME);
		b.append(name);
		b.append(COMMENT_SEP);
		b.append(Isoform.GID);
		b.append(name);
		b.append(GENNAME);
		b.append(name);
		b.append(Isoform.ID);
		b.append(name);
		
		return b.toString();
	}
	
	/**
	 * returns the name for the gff file of multiple exons
	 * @param exon
	 * @return
	 */
	protected String getName(Exon ...exon) {
		StringBuffer name = new StringBuffer();
		for (int i = 0; i < exon.length; i++) {
			name.append(exon[i].toString());
			if(i + 1 < exon.length) {
				name.append(SEP_NAME);
			}
		}
		return name.toString();
	}
	
	/**
	 * Writes the isoforms to a file or the buffer
	 * @param geneHeader
	 * @param name
	 * @param isoform
	 * @throws IOException 
	 */
	protected void write(String geneHeader, String name, Isoform ...isoform) throws IOException {
		// do nothing in case of silent mode because only the name should be fetched
		if(!this.SILENT_MODE) {
			int ii = 1;
			
			if(this.BW != null) {
				this.BW.write(geneHeader);
				this.BW.write(NEWLINE);
				
				for(Isoform i : isoform) {
					this.BW.write(i.toGFF3(name.toString(), ii));
					this.BW.write(NEWLINE);
					ii++;
				}
			}
			else {
				
				this.BW.append(geneHeader);
				this.BW.append(NEWLINE);
				
				for(Isoform i : isoform) {
					this.BW.append(i.toGFF3(name.toString(), ii));
					this.BW.append(NEWLINE);
					ii++;
				}
			}
		}
	}
	
	/**
	 * checks, if a event is significant or not based on the FDR value
	 * @param split
	 * @return
	 */
	protected boolean isSignificant(String split[]) {
		if(this.SILENT_MODE && this.FDR != null) {
			return Double.parseDouble(split[this.NAME2ID.get(NAME_FDR)]) <= this.FDR;
		}
		return false;
	}
	
	/**
	 * significant event names or emtpy list if no valid FDR was set
	 * @return
	 */
	public ArrayList<String> getSigEvents() {
		return new ArrayList<String>(this.EVENT_NAMES);
	}
	
	/****************************** ABSTRACT METHODS *****************************/
	
	/**
	 * Processes a single line of the file and adds the gff3 annotation to the StringBuffer
	 * @param split
	 * @return name of the event
	 * @throws IOException 
	 */
	protected abstract String getEvent(String split[]) throws IOException;
	
	/**
	 * Returns the type of the event
	 * @return
	 */
	public abstract String getType();
}
