package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Parses the SE Annotation
 * @author Michael Kluge
 *
 */
public class SEParser extends MatsOutputParser {
	
	private static final String TYPE = "SE";
	
	private static final String GENE_NAME = "GeneID";
	private static final String CHR_NAME = "chr";
	private static final String STRAND_NAME = "strand";
	
	private static final String EXON_SE_START_NAME = "exonStart_0base";
	private static final String EXON_SE_END_NAME = "exonEnd";
	private static final String EXON_UP_START_NAME = "upstreamES";
	private static final String EXON_UP_END_NAME = "upstreamEE";
	private static final String EXON_DOWN_START_NAME = "downstreamES";
	private static final String EXON_DOWN_END_NAME = "downstreamEE";
	
	static {
		ArrayList<String> names = new ArrayList<String>();
		NAMES.put(TYPE, names);
		
		names.add(GENE_NAME);
		names.add(CHR_NAME);
		names.add(STRAND_NAME);
		names.add(EXON_SE_START_NAME);
		names.add(EXON_SE_END_NAME);
		names.add(EXON_UP_START_NAME);
		names.add(EXON_UP_END_NAME);
		names.add(EXON_DOWN_START_NAME);
		names.add(EXON_DOWN_END_NAME);
	}
		

	public SEParser(String path, BufferedWriter bw) throws IOException {
		super(path, bw);
	}

	@Override
	protected void getEvent(String split[]) throws IOException {
		
		String gene = split[this.NAME2ID.get(GENE_NAME)];
		String chr = split[this.NAME2ID.get(CHR_NAME)];
		boolean plusStrand = (split[this.NAME2ID.get(STRAND_NAME)]).equals(Exon.PLUS_STRAND);
		int SEs = Integer.parseInt(split[this.NAME2ID.get(EXON_SE_START_NAME)]);
		int SEe = Integer.parseInt(split[this.NAME2ID.get(EXON_SE_END_NAME)]);
		int UPs = Integer.parseInt(split[this.NAME2ID.get(EXON_UP_START_NAME)]);
		int UPe = Integer.parseInt(split[this.NAME2ID.get(EXON_UP_END_NAME)]);
		int DOs = Integer.parseInt(split[this.NAME2ID.get(EXON_DOWN_START_NAME)]);
		int DOe = Integer.parseInt(split[this.NAME2ID.get(EXON_DOWN_END_NAME)]);
		
		Exon up = new Exon(chr, UPs, UPe, plusStrand);
		Exon se = new Exon(chr, SEs, SEe, plusStrand);
		Exon down = new Exon(chr, DOs, DOe, plusStrand);
		
		Isoform i1 = new Isoform(this.getType(), up, se, down);
		Isoform i2 = new Isoform(this.getType(), up, down);
		
		String name = this.getName(up, se, down);
		String geneLine = this.getGeneLine(gene, chr, up, se, down);
		
		this.write(geneLine, name, i1, i2);
	}

	@Override
	protected String getType() {
		return TYPE;
	}
}