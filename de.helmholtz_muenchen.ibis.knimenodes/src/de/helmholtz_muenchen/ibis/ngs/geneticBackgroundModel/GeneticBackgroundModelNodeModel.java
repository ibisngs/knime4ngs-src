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


package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.Identifier.EntityIdentifier;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.MatrixSummary;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.BioEntity;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;
import de.helmholtz_muenchen.ibis.utils.ngs.VEPAnnotationParser;

/**
 * This is the model implementation of GeneticBackgroundModel.
 * 
 *
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeModel extends NodeModel {
	
	//configuration keys
	static final String CFGKEY_GTF_AFF = "gtf_aff";
	static final String CFGKEY_AC = "ac";
	static final String CFGKEY_AN = "an";
	static final String CFGKEY_RESOLUTION = "resolution";
	static final String [] RESOLUTION = {"gene_id", "gene_symbol", "transcript_id"};
	
	//settings models
	private final SettingsModelString m_ac = new SettingsModelString(CFGKEY_AC,"AC");
	private final SettingsModelString m_an = new SettingsModelString(CFGKEY_AN,"AN");
	private final SettingsModelString m_resolution = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_RESOLUTION, GeneticBackgroundModelNodeModel.RESOLUTION[0]);
	
	private boolean use_genotypes;
	
	//output col names
	public static final String OUT_COL1 = "Path2BackgroundModel";
	
	NodeLogger LOGGER = NodeLogger.getLogger(GeneticBackgroundModelNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected GeneticBackgroundModelNodeModel() {
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//get and check input file
    	String infile = inData[0].iterator().next().getCell(0).toString();
    	
    	if(Files.notExists(Paths.get(infile))) {
    		throw new InvalidSettingsException("Input file does not exist!");
    	}

    	//check resolution
    	String ending = m_resolution.getStringValue();
    	BioEntity e;
    	if(ending.equals("gene_id")) e = BioEntity.GENE_ID;
    	else if(ending.equals("gene_symbol")) e = BioEntity.GENE_SYMBOL;
    	else e = BioEntity.TRANSCRIPT_ID;

    	//compute frequencies according to input file type (VCF or matrix)
    	HashMap<String, Double> gene_frequency;
    	if(use_genotypes) {
    		
    		if(!infile.endsWith("matrix.tsv")) {
    			throw new InvalidSettingsException("The input file is neither a VCF nor a matrix file!");
    		}
    		
    		MatrixSummary ms = new MatrixSummary(infile);
    		if(e == BioEntity.TRANSCRIPT_ID) {
    			gene_frequency = ms.getFrequencies();
    		} else {
    			gene_frequency = ms.getFrequencies(new EntityIdentifier(e));
    		}
    		ending += ".model_gtf.tsv";
    	} else {
    		VCFFile vcf_it = new VCFFile(infile);
        	String vep_header = vcf_it.getInfoHeader(VEPAnnotationParser.ANN_ID);
        	
        	if(vep_header == null) {
        		throw new InvalidSettingsException("No VEP annotations found!");
        	}
        	
        	AnnotationParser parser = new VEPAnnotationParser(vep_header);
        	gene_frequency = fillAF(vcf_it, m_ac.getStringValue(), m_an.getStringValue(), parser, e);
    		ending += ".model_aff.tsv";
    	}
    	
		String outfile = IO.replaceFileExtension(infile, ending);
    	writeModel(outfile, gene_frequency);
		
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
        return new BufferedDataTable[]{outTable};
    }

	private HashMap<String, Double> fillAF(VCFFile vcf_it, String ac_id, String an_id, AnnotationParser parser, BioEntity ent) {
		HashMap<String, Double> result = new HashMap<>();
		VCFVariant var;
		String ac, an, csq;
		while(vcf_it.hasNext()) {
			var = vcf_it.next();
			ac = var.getInfoField(ac_id);
			an = var.getInfoField(an_id);
			csq = var.getInfoField(parser.getAnnId());
			
			if(ac == null || an == null) {
				LOGGER.error("The INFO fields "+ac_id+ " and "+an_id+ " have not been found for variant on chr "+var.getChrom()+" at position "+var.getPos()+"!");
				continue;
			}
			
			if(csq == null) {
				LOGGER.error("No annotations have been found for variant on chr "+var.getChrom()+" at position "+var.getPos()+"!");
				continue;
			}
			
			HashMap<String, HashSet<Integer>> gene2allele_num = parser.getEntity2AlleleIds(csq, ent);
			String [] acs = ac.split(",");
			
			//compute frequency of being unaffected for each gene 
			double ac_abs, af, not_aff;
			
			for(String g: gene2allele_num.keySet()) {
				ac_abs = 0;
				for(int i: gene2allele_num.get(g)) {
					ac_abs += Integer.parseInt(acs[i-1]);
				}
				
				af = (double)ac_abs/Double.parseDouble(an);
				not_aff = Math.pow(1.0-af,2);
				if(result.containsKey(g)) {
					not_aff = result.get(g) * not_aff;
				}
				result.put(g, not_aff);
			}
			
		}
		
		for(String g: result.keySet()) {
			result.put(g, 1.0 - result.get(g));
		}
		
		return result;
	}
	
    private void writeModel(String outfile, HashMap<String, Double> gene_frequency) {
		
    	Charset c = Charset.forName("UTF-8");
    	BufferedWriter bw;
    	try {
			bw = Files.newBufferedWriter(Paths.get(outfile), c);
			bw.write(m_resolution.getStringValue()+"\tvariant_freq");
			bw.newLine();
			for(String s: gene_frequency.keySet()) {
				double freq = gene_frequency.get(s);
				if(freq > 0.0) {
					bw.write(s+"\t"+gene_frequency.get(s));
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {


    	if(inSpecs[0].getColumnSpec(0).getType().toString().equals("VCFCell")) {
    			use_genotypes = false;
    	} else {
    		use_genotypes = true;
    	}
    	
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_ac.saveSettingsTo(settings);
    	m_an.saveSettingsTo(settings);
    	m_resolution.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_ac.loadSettingsFrom(settings);
    	m_an.loadSettingsFrom(settings);
    	m_resolution.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_ac.validateSettings(settings);
    	m_an.validateSettings(settings);
    	m_resolution.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    }

}

