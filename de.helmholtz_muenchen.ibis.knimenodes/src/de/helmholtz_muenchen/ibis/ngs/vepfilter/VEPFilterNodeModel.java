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
package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
/**
 * This is the model implementation of VEPFilter.
 * 
 *
 * @author Tim Jeske
 */
public class VEPFilterNodeModel extends HTExecutorNodeModel {
   
	static final String DEFAULT_TERM = "NONE";
	
	static final String [] SO_TERMS = {
			DEFAULT_TERM,
			"transcript_ablation",
			"splice_acceptor_variant",
			"splice_donor_variant",
	    	"stop_gained",
	    	"frameshift_variant",
	    	"stop_lost",
	    	"start_lost",
	    	"transcript_amplification",
	    	"inframe_insertion",
	    	"inframe_deletion",
	    	"missense_variant",
	    	"protein_altering_variant",
	    	"splice_region_variant",
	    	"incomplete_terminal_codon_variant",
	    	"stop_retained_variant",
	    	"synonymous_variant",
	    	"coding_sequence_variant",
	    	"mature_miRNA_variant",
	    	"5_prime_UTR_variant",
	    	"3_prime_UTR_variant",
	    	"non_coding_transcript_exon_variant",
	    	"intron_variant",
	    	"NMD_transcript_variant",
	    	"non_coding_transcript_variant",
	    	"upstream_gene_variant",
	    	"downstream_gene_variant",
	    	"TFBS_ablation",
	    	"TFBS_amplification",
	    	"TF_binding_site_variant",
	    	"regulatory_region_ablation",
	    	"regulatory_region_amplification",
	    	"feature_elongation",
	    	"regulatory_region_variant",
	    	"feature_truncation",
	    	"intergenic_variant"
	};
	
	static final String [] LOF_TERMS = {SO_TERMS[1], SO_TERMS[2], SO_TERMS[3], SO_TERMS[4]};
	
	//configuration keys
    static final String CFGKEY_VEP_SCRIPT = "vepscript";
    static final String CFGKEY_SO_TERM = "so_term"; 
    static final String CFGKEY_TERM_LIST = "term_list";
    static final String CFGKEY_FILTER = "filter";    
    static final String CFGKEY_OUTFOLDER = "outfolder";
	static final String CFGKEY_OVERWRITE = "overwrite";
	
	//settings models
	private final SettingsModelString m_vepscript = new SettingsModelString(CFGKEY_VEP_SCRIPT,"-");
    private final SettingsModelOptionalString m_filter = new SettingsModelOptionalString(CFGKEY_FILTER,"",false);
	private final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
	private final SettingsModelBoolean m_overwrite = new SettingsModelBoolean(CFGKEY_OVERWRITE, false);
	private final SettingsModelStringArray m_chosen_terms = new SettingsModelStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST, new String[]{VEPFilterNodeModel.DEFAULT_TERM});
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VEPFilterNodeModel.class);
	
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
	private int vcf_index;
	private String vep_script;
	
    protected VEPFilterNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1), 1);
    	
    	boolean use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
    	
    	if(use_hte) {
    		m_overwrite.setBooleanValue(true);
    		m_overwrite.setEnabled(false);
    	}
    	
    	addSetting(m_filter);
    	addSetting(m_overwrite);
    	addSetting(m_outfolder);
    	addSetting(m_vepscript);
    	addSetting(m_chosen_terms);
    	
    	addPrefPageSetting(m_vepscript, IBISKNIMENodesPlugin.VEP_FILTER);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	String infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(infile))) {
    		throw new InvalidSettingsException("Input VCF does not exist!");
    	}
    	
    	String outfile = filterAnnotations(infile, exec);
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	return new BufferedDataTable[]{outTable};
    }
    
    private String filterAnnotations(String infile, ExecutionContext exec) throws Exception {
    	
		String outfile;
		if (infile.endsWith(".gz")) {
			outfile = infile.replace(".vcf.gz", ".VEPfiltered.vcf");
		} else {
			outfile = IO.replaceFileExtension(infile,".VEPfiltered.vcf");
		}
		
		String outfolder = IO.processFilePath(m_outfolder.getStringValue());
    	if(outfolder!=null && !outfolder.equals("")) {
    		outfile = outfolder + System.getProperty("file.separator") + new File(outfile).getName();
    	}

		LOGGER.debug("CHOSEN TERMS: " + m_chosen_terms.getStringArrayValue());
		LOGGER.debug("FILTER: " + m_filter.getStringValue());

		LOGGER.info("Prepare command for filter_vep.pl");
		ArrayList<String> cmd = new ArrayList<>();
		cmd.add("perl");
		cmd.add(vep_script);
		cmd.add("-i");
		cmd.add(infile);
		cmd.add("-o");
		cmd.add(outfile);

		String filter_terms = "";
		String[] terms = m_chosen_terms.getStringArrayValue();
		if (terms.length > 0 && !terms[0].equals(DEFAULT_TERM)) {
			cmd.add("--filter");
			filter_terms += "Consequence is " + terms[0];
			for (int i = 1; i < terms.length; i++) {
				filter_terms += " or Consequence is " + terms[i];
			}
			cmd.add(filter_terms);
		}

		if(m_filter.isActive()) {
			String[] filters = m_filter.getStringValue().split(",");
			if (filters.length > 0) {
				for (String f : filters) {
					if (f.equals(""))
						continue;
					cmd.add("--filter");
					f = f.trim();
					cmd.add(f);
				}
			}
		}

		cmd.add("--only_matched");
		
		if(m_overwrite.getBooleanValue()) {
			cmd.add("--force_overwrite");
		}

		String[] cmd_array = new String[cmd.size()];
		for (int i = 0; i < cmd.size(); i++) {
			cmd_array[i] = cmd.get(i);
		}

		File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);
		
		String perl5lib_variable = "PERL5LIB="+System.getenv("PERL5LIB");
		
		super.executeCommand(cmd_array, outfile, exec, new String[]{perl5lib_variable}, lockFile, null, null, null ,null, null);
		
		return outfile;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	super.updatePrefs();
    	vep_script = IO.processFilePath(m_vepscript.getStringValue());
    	
    	if(CompatibilityChecker.inputFileNotOk(vep_script, false)) {
    		throw new InvalidSettingsException("Invalid path to filter_vep.pl script!");
    	}
    	
    	if(VEPFilterNodeDialog.getUseMainInputColBool()){
    		vcf_index = inSpecs[0].findColumnIndex(VEPFilterNodeDialog.getMainInputCol1());
    		if(!inSpecs[0].getColumnSpec(vcf_index).getType().toString().equals("VCFCell")){
    			vcf_index = -1;
    		}
    	} else {
    		vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
    	}
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the output table!");
    	}
    	
    	String [] terms = m_chosen_terms.getStringArrayValue();
    	for(String s: terms) {
    		if(s.equals(DEFAULT_TERM) && terms.length > 1) {
    			throw new InvalidSettingsException("SO terms including "+DEFAULT_TERM + " are selected!");
    		}
    	}
    	
    	if(terms[0].equals(DEFAULT_TERM) && (!m_filter.isActive() || m_filter.getStringValue().equals(""))) {
    		throw new InvalidSettingsException("No valid filters given!");
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }
}

