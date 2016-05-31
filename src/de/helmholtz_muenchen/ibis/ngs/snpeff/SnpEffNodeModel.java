package de.helmholtz_muenchen.ibis.ngs.snpeff;

import java.io.File;
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
	private final SettingsModelString m_database = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATABASE,"");
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

	
	//The Output Col Names
	public static final String OUT_COL1 = "outputVCF";
	
	private int vcf_index;
	
	/**
     * Constructor for the node model.
     */
    protected SnpEffNodeModel() {
    
        super(OptionalPorts.createOPOs(1),OptionalPorts.createOPOs(1));
        
    	addSetting(m_snpeff_bin);
    	addSetting(m_database);
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
    			(FileCell) FileCellFactory.create(out_file)};
    	
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

    	vcf_index = -1;
    	for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
    	
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
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }
}