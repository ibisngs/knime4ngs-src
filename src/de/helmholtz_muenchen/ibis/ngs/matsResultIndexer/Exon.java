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

/**
 * Stores information about exons
 * @author Michael Kluge
 *
 */
public class Exon {
	
	private static final String SEP = ":";
	public static final char PLUS_STRAND = '+';
	public static final char NEG_STRAND = '-';
	
	private final String CHR;
	private final int START;
	private final int END;
	private char STRAND;
	
	/**
	 * Constructor
	 * @param chr
	 * @param start
	 * @param end
	 * @param plusStrand
	 */
	public Exon(String chr, int start, int end, boolean plusStrand) {
		this.CHR = chr;
		this.START = start;
		this.END = end;
		
		if(plusStrand) this.STRAND = PLUS_STRAND;
		else this.STRAND = NEG_STRAND;
	}

	public char getStrand() {
		return STRAND;
	}

	public String getChr() {
		return CHR;
	}

	public int getStart() {
		return START;
	}

	public int getEnd() {
		return END;
	}
	
	public String toString() {
		return CHR + SEP + START + SEP + END + SEP + STRAND; 
	}
}
