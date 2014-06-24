package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Parses the A5SS Annotation
 * @author Michael Kluge
 *
 */
public class A5SSParser extends MatsOutputParser {
	
	private static final String TYPE = "A5SS";
	
	private static final String GENE_NAME = "GeneID";
	private static final String CHR_NAME = "chr";
	private static final String STRAND_NAME = "strand";
	
	private static final String EXON_LONG_START_NAME = "longExonStart_0base";
	private static final String EXON_LONG_END_NAME = "longExonEnd";
	private static final String EXON_SHORT_START_NAME = "shortES";
	private static final String EXON_SHORT_END_NAME = "shortEE";
	private static final String EXON_FLANKING_START_NAME = "flankingES";
	private static final String EXON_FLANKING_END_NAME = "flankingEE";
	
	static {
		ArrayList<String> names = new ArrayList<String>();
		NAMES.put(TYPE, names);
		
		names.add(GENE_NAME);
		names.add(CHR_NAME);
		names.add(STRAND_NAME);
		names.add(EXON_LONG_START_NAME);
		names.add(EXON_LONG_END_NAME);
		names.add(EXON_SHORT_START_NAME);
		names.add(EXON_SHORT_END_NAME);
		names.add(EXON_FLANKING_START_NAME);
		names.add(EXON_FLANKING_END_NAME);
	}
		

	public A5SSParser(String path, BufferedWriter bw) throws IOException {
		super(path, bw);
	}

	@Override
	protected void getEvent(String split[]) throws IOException {
		
		String gene = split[this.NAME2ID.get(GENE_NAME)];
		String chr = split[this.NAME2ID.get(CHR_NAME)];
		boolean plusStrand = (split[this.NAME2ID.get(STRAND_NAME)]).equals(Exon.PLUS_STRAND);
		int LongEs = Integer.parseInt(split[this.NAME2ID.get(EXON_LONG_START_NAME)]);
		int LongEe = Integer.parseInt(split[this.NAME2ID.get(EXON_LONG_END_NAME)]);
		int ShortEs = Integer.parseInt(split[this.NAME2ID.get(EXON_SHORT_START_NAME)]);
		int ShortEe = Integer.parseInt(split[this.NAME2ID.get(EXON_SHORT_END_NAME)]);
		int FlankEs = Integer.parseInt(split[this.NAME2ID.get(EXON_FLANKING_START_NAME)]);
		int FlankEe = Integer.parseInt(split[this.NAME2ID.get(EXON_FLANKING_END_NAME)]);
		
		Exon longE = new Exon(chr, LongEs, LongEe, plusStrand);
		Exon shortE = new Exon(chr, ShortEs, ShortEe, plusStrand);
		Exon flankE = new Exon(chr, FlankEs, FlankEe, plusStrand);
		
		Isoform i1 = new Isoform(this.getType(), longE, flankE);
		Isoform i2 = new Isoform(this.getType(), shortE, flankE);
		
		String name = this.getName(longE, shortE, flankE);
		String geneLine = this.getGeneLine(gene, chr, longE, shortE, flankE);
		
		this.write(geneLine, name, i1, i2);
	}

	@Override
	protected String getType() {
		return TYPE;
	}
}