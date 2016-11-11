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
package de.helmholtz_muenchen.ibis.ngs.snpeff;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of SnpEff.
 * 
 *
 * @author Sebastian Kopetzky
 */
public class SnpEffNodeModel extends HTExecutorNodeModel {
	
	static final String CFGKEY_SNPEFF_BIN = "snpeff_binary";
	static final String CFGKEY_DATABASE = "database";
	static final String CFGKEY_USE_CONFIG = "use_config";
	static final String CFGKEY_DATA_DIR = "database_dir";
	static final String CFGKEY_USEBEDFILE="usebedfile";
	static final String CFGKEY_BED_FILE = "bed_file";
	static final String CFGKEY_MEM = "memory";
	static final String CFGKEY_OPT_FLAGS ="opt_flags";
	
	//results filter
	static final String CFGKEY_NO_DOWNSTREAM = "no_downstream";
	static final String CFGKEY_NO_INTERGENIC = "no_intergenic";
	static final String CFGKEY_NO_INTRONIC = "no_intronic";
	static final String CFGKEY_NO_UPSTREAM = "no_upstream";
	static final String CFGKEY_NO_UTR = "no_utr";

	private final SettingsModelString m_snpeff_bin = new SettingsModelString(SnpEffNodeModel.CFGKEY_SNPEFF_BIN,"");
	private final SettingsModelBoolean m_use_config = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_USE_CONFIG, true);
	private final SettingsModelString m_database = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATABASE,"");
	private final SettingsModelString m_database_dir = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATA_DIR, "");
	private final SettingsModelBoolean m_usebedfile = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_USEBEDFILE, false);
	private final SettingsModelString m_bed_file = new SettingsModelString(SnpEffNodeModel.CFGKEY_BED_FILE,"");
	private final SettingsModelIntegerBounded m_memory = new SettingsModelIntegerBounded(SnpEffNodeModel.CFGKEY_MEM, 4, 1, Integer.MAX_VALUE);
    private final SettingsModelOptionalString m_opt_flags = new SettingsModelOptionalString(SnpEffNodeModel.CFGKEY_OPT_FLAGS,"",false);
	
	//results filter
	private final SettingsModelBoolean m_no_downstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_DOWNSTREAM, false);
	private final SettingsModelBoolean m_no_intergenic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTERGENIC, false);
	private final SettingsModelBoolean m_no_intronic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTRONIC, false);
	private final SettingsModelBoolean m_no_upstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UPSTREAM, false);
	private final SettingsModelBoolean m_no_utr = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UTR, false);

	private static final NodeLogger LOGGER = NodeLogger.getLogger(SnpEffNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "outputVCF";
	
	private int vcf_index;
	
	/**
     * Constructor for the node model.
     */
    protected SnpEffNodeModel() {
    
        super(OptionalPorts.createOPOs(1),OptionalPorts.createOPOs(1));
        
    	addSetting(m_snpeff_bin);
    	addSetting(m_use_config);
    	addSetting(m_database);
    	addSetting(m_database_dir);
    	addSetting(m_usebedfile);
    	addSetting(m_bed_file);
    	addSetting(m_memory);
    	addSetting(m_opt_flags);
    	
    	//results filter
    	addSetting(m_no_downstream);
    	addSetting(m_no_intergenic);
    	addSetting(m_no_intronic);
    	addSetting(m_no_upstream);
    	addSetting(m_no_utr);
        
        m_bed_file.setEnabled(false);
        m_database_dir.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	/*Get parameters*/
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(CompatibilityChecker.inputFileNotOk(vcf_infile)) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String snpEffBin = m_snpeff_bin.getStringValue();
    	
    	//get database dir 
    	String db_dir;
    	if(m_use_config.getBooleanValue()) {
    		db_dir = getDataDirFromConfig(getConfigFile());
    	} else {
    		db_dir = m_database_dir.getStringValue();
    	}
    	LOGGER.debug("Database dir: " +db_dir);
    	
    	//check if database exists and download it
    	String database = "";
    	if(!db_dir.endsWith(File.separator)) {
    		database = db_dir + File.separator;
    	}
    	database += m_database.getStringValue();
    	LOGGER.debug("Database path: "+ database);
    	
    	if(!Files.exists(Paths.get(database))) {
    		ArrayList<String> command = new ArrayList<String>();
    		command.add("java");
    		command.add("-jar " + snpEffBin + " download");
    		command.add(m_database.getStringValue());
        	
        	LOGGER.debug("Download database "+m_database.getStringValue());
        	String lockFile = new File(vcf_infile).getParent() + File.separator + "snpEff_download" + SuccessfulRunChecker.LOCK_ENDING;
        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));
        	LOGGER.debug("Database download finished!");
    	}
    	
    	String out_file = IO.replaceFileExtension(vcf_infile, "snpEff.vcf");
    	String stats_file = IO.replaceFileExtension(vcf_infile, "summary.html");
    	   	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-Xmx"+m_memory.getIntValue()+"G -jar "+snpEffBin);
        command.add("-s "+stats_file);
    	
    	//Result filter options
    	if(m_usebedfile.getBooleanValue()){
    		command.add("-fi " + m_bed_file.getStringValue());
    	}
    	if(m_no_downstream.getBooleanValue()){
    		command.add("-no-downstream");
    	}
    	if(m_no_intergenic.getBooleanValue()){
    		command.add("-no-intergenic");
    	}
    	if(m_no_intronic.getBooleanValue()){
    		command.add("-no-intron");
    	}
    	if(m_no_upstream.getBooleanValue()){
    		command.add("-no-upstream");
    	}
    	if(m_no_utr.getBooleanValue()){
    		command.add("-no-utr");
    	}
    	
    	command.add(m_opt_flags.getStringValue());
    	command.add("-v "+m_database.getStringValue());
		command.add("-dataDir "+ db_dir);
    	command.add(vcf_infile);
    	
    	/**Execute**/
    	String lockFile = out_file + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile),out_file);
    	
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(out_file)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_snpeff_bin.getStringValue())) {
    		throw new InvalidSettingsException("Set a valid path to the snpEff directory!");
    	}
    	
    	String db = m_database.getStringValue();
    	if(db.equals("") || db == null) {
    		throw new InvalidSettingsException("Specify an annotation database!");
    	}
    	
    	if(m_use_config.getBooleanValue() && CompatibilityChecker.inputFileNotOk(getConfigFile())) {
    		throw new InvalidSettingsException("Config file not found or invalid!");
    	}
    	
    	if(!m_use_config.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_database_dir.getStringValue(),false)) {
    		throw new InvalidSettingsException("Specify a valid path to the database directory!");
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }
    
    private String getConfigFile() {
    	return new File(m_snpeff_bin.getStringValue()).getParent() + File.separator + "snpEff.config";
    }
    
    private String getDataDirFromConfig(String configPath) throws InvalidSettingsException {
    	
    	String db_dir = "";
    	BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(configPath));
			String line;
	    	while((line = br.readLine())!= null) {
	    		if(line.startsWith("data.dir")) {
	    			db_dir = line.split("=")[1].trim();
	    			if(db_dir.startsWith(".")) {
	    				db_dir = new File(configPath).getParent() + File.separator + db_dir.split(File.separator)[1];
	    			}
	    			br.close();
	    			return db_dir;
	    		}
	    	}
	    	br.close();
		} catch (IOException e) {
			throw new InvalidSettingsException("Config file not in the given snpEff directory!");
		} 
    	
    	return db_dir;
    }
}