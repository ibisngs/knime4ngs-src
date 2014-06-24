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
