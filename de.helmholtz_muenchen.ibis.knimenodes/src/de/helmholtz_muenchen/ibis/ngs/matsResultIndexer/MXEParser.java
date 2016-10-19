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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Parses the MXE Annotation
 * @author Michael Kluge
 *
 */
public class MXEParser extends MatsOutputParser {
	
	private static final String TYPE = "MXE";
	
	private static final String GENE_NAME = "GeneID";
	private static final String CHR_NAME = "chr";
	private static final String STRAND_NAME = "strand";
	
	private static final String EXON_FIRST_START_NAME = "1stExonStart_0base";
	private static final String EXON_FIRST_END_NAME = "1stExonEnd";
	private static final String EXON_SECOND_START_NAME = "2ndExonStart_0base";
	private static final String EXON_SECOND_END_NAME = "2ndExonEnd";
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
		names.add(EXON_FIRST_START_NAME);
		names.add(EXON_FIRST_END_NAME);
		names.add(EXON_SECOND_START_NAME);
		names.add(EXON_SECOND_END_NAME);
		names.add(EXON_UP_START_NAME);
		names.add(EXON_UP_END_NAME);
		names.add(EXON_DOWN_START_NAME);
		names.add(EXON_DOWN_END_NAME);
	}
	
	public MXEParser(String path, BufferedWriter bw, Double FDR) throws IOException {
		super(path, bw, FDR);
	}

	@Override
	protected String getEvent(String split[]) throws IOException {
		
		String gene = split[this.NAME2ID.get(GENE_NAME)];
		String chr = split[this.NAME2ID.get(CHR_NAME)];
		boolean plusStrand = (split[this.NAME2ID.get(STRAND_NAME)]).toCharArray()[0] == Exon.PLUS_STRAND;
		int firstExonS = Integer.parseInt(split[this.NAME2ID.get(EXON_FIRST_START_NAME)]);
		int firstExonE = Integer.parseInt(split[this.NAME2ID.get(EXON_FIRST_END_NAME)]);
		int secondExonS = Integer.parseInt(split[this.NAME2ID.get(EXON_SECOND_START_NAME)]);
		int secondExonE = Integer.parseInt(split[this.NAME2ID.get(EXON_SECOND_END_NAME)]);
		int UPs = Integer.parseInt(split[this.NAME2ID.get(EXON_UP_START_NAME)]);
		int UPe = Integer.parseInt(split[this.NAME2ID.get(EXON_UP_END_NAME)]);
		int DOs = Integer.parseInt(split[this.NAME2ID.get(EXON_DOWN_START_NAME)]);
		int DOe = Integer.parseInt(split[this.NAME2ID.get(EXON_DOWN_END_NAME)]);
		
		Exon up = new Exon(chr, UPs, UPe, plusStrand);
		Exon first = new Exon(chr, firstExonS, firstExonE, plusStrand);
		Exon second = new Exon(chr, secondExonS, secondExonE, plusStrand);
		Exon down = new Exon(chr, DOs, DOe, plusStrand);
		
		Isoform i1 = new Isoform(this.getType(), up, first, down);
		Isoform i2 = new Isoform(this.getType(), up, second, down);
		
		String name = this.getName(up, first, second, down);
		String geneLine = this.getGeneLine(gene, chr, up, first, second, down);
		
		this.write(geneLine, name, i1, i2);
		
		if(this.isSignificant(split)) {
			return name;
		}
		return null;
	}

	@Override
	public String getType() {
		return TYPE;
	}
}