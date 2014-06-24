package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Stores the exons of a isoform and is able to output the annotation in a gff format which is 
 * accepted by MISO.
 * @author Michael Kluge
 *
 */
public class Isoform {
	
	private static final int ASCII_OFFSET = 64;
	public static final String GID = "gid=";
	public static final String NAME_SEP = ".";
	public static final String PARENT = ";Parent=";
	public static final String ID = ";ID=";
	public static final String MRNA = "mRNA";
	public static final String EXON = "exon";
	public static String TAB = "\t";
	public static String NEWLINE = "\n";
	public static String SPACER = ".";
	private final ArrayList<Exon> EXONS = new ArrayList<Exon>();
	private final String CHR;
	private final String TYPE;
	
	/**
	 * Constructor
	 * @param type
	 * @param exon
	 */
	public Isoform(String type, Exon ... exon) {
		this.TYPE = type;
		// add exons
		for (int i = 0; i < exon.length; i++) {
			this.EXONS.add(exon[i]);
		}
		this.CHR = this.getFirst().getChr();
	}
	
	/**
	 * Ouputs the isoform in gff format which can be read by MISO
	 * @param geneName identifier of the gene
	 * @param isoformNumber number of the isoform of the gene
	 * @return
	 */
	public String toGFF3(String geneName, int isoformNumber) {
		StringBuffer b = new StringBuffer();
		
		// mRNA line
		b.append(this.CHR);
		b.append(TAB);
		b.append(this.TYPE);
		b.append(TAB);
		b.append(MRNA);
		b.append(TAB);
		b.append(this.getFirst().getStart());
		b.append(TAB);
		b.append(this.getLast().getEnd());
		b.append(TAB);
		b.append(SPACER);
		b.append(TAB);
		b.append(this.getFirst().getStrand());
		b.append(TAB);
		b.append(SPACER);
		b.append(TAB);
		// add comment stuff;
		b.append(GID);
		b.append(geneName);
		b.append(ID);
		b.append(geneName);
		b.append(NAME_SEP);
		b.append((char) (ASCII_OFFSET+isoformNumber));
		b.append(PARENT);
		b.append(geneName);
		b.append(NEWLINE);
		
		// exon line
		int i=1;
		for(Iterator<Exon> it = this.EXONS.iterator(); it.hasNext(); ) {
			Exon e = it.next();

			b.append(this.CHR);
			b.append(TAB);
			b.append(this.TYPE);
			b.append(TAB);
			b.append(EXON);
			b.append(TAB);
			b.append(e.getStart());
			b.append(TAB);
			b.append(e.getEnd());
			b.append(TAB);
			b.append(SPACER);
			b.append(TAB);
			b.append(this.getFirst().getStrand());
			b.append(TAB);
			b.append(SPACER);
			b.append(TAB);
			// add comment stuff;
			b.append(GID);
			b.append(geneName);
			b.append(ID);
			b.append(geneName);
			b.append(NAME_SEP);
			b.append((char) (ASCII_OFFSET + isoformNumber));
			b.append(NAME_SEP);
			b.append(i);
			b.append(PARENT);
			b.append(geneName);
			b.append(NAME_SEP);
			b.append((char) (ASCII_OFFSET + isoformNumber));
			
			if(it.hasNext()) b.append(NEWLINE);
			i++;
		}
		return b.toString();
	}
	
	/**
	 * returns the first exon of the isoform
	 * @return
	 */
	public Exon getFirst() {
		return this.EXONS.get(0);
	}
	
	/**
	 * returns the last exon of the isoform
	 * @return
	 */
	public Exon getLast() {
		return this.EXONS.get(EXONS.size() - 1);
	}
}
