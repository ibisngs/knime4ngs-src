package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;



/**
 * This is the model implementation of DepthOfCoverage.
 * 
 *
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class DepthOfCoverageNodeModel extends GATKNodeModel {
    
	 protected static final NodeLogger logger = NodeLogger.getLogger(DepthOfCoverageNodeModel.class);
	
	/**
	 * Config Keys
	 */
	
	public static final String CFGKEY_PATH2BED = "bed";
	public static final String CFGKEY_BED_FILE_CHECKBOX = "BED_FILE_CHECKBOX";
	public static final String CFGKEY_INFILE = "infile";
	public static final String CFGKEY_EXTRAFILTERS = "extrafilters";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	/**
	 * Node Models
	 */
	private final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
    private final SettingsModelBoolean m_bed_file_checkbox = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
	private final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
	private final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"DoC");
	private final SettingsModelString m_infile = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_INFILE,"");
	
	private String OUTFILE, LOCKFILE;

    /**
     * Constructor for the node model.
     */
    protected DepthOfCoverageNodeModel() {
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    @Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
    	
    	String infile = "";
    	ArrayList<String> command = new ArrayList<String>();
    	
		if (inData[0]!=null) {
			infile = inData[0].iterator().next().getCell(0).toString();
		} else {
			infile = m_infile.getStringValue();
		}
		
		if(infile.equals("") || Files.notExists(Paths.get(infile))) {
			throw new InvalidSettingsException("No infile specified!");
		}
		
		String fileSuffix = m_filesuffix.getStringValue();
		if(fileSuffix.equals("")) {
			fileSuffix = "DoC";
		}
		
		this.OUTFILE = IO.replaceFileExtension(infile, fileSuffix);
		this.LOCKFILE = IO.replaceFileExtension(infile, SuccessfulRunChecker.LOCK_ENDING);
		command.add("-I "+infile);

		if(m_bed_file_checkbox.getBooleanValue()){
			if(m_path2bed.getStringValue().equals("") || Files.notExists(Paths.get(m_path2bed.getStringValue()))) {
				setWarningMessage("Specify valid bed file!");
			} else {
				command.add("-L "+m_path2bed.getStringValue());
			}
    	}
		
		if(!m_extrafilters.getStringValue().equals("")){
        	String[] filters = m_extrafilters.getStringValue().split(",");
        	for(String filter : filters){
        		command.add(" --read_filter "+filter);
        	}
    	}
		
		DepthOfCoverageNodeModel.logger.info("Running GATK DepthOfCoverage...");
    	DepthOfCoverageNodeModel.logger.info("Log files can be found in "+OUTFILE+".out.log and "+OUTFILE+".err.log");
		
		return StringUtils.join(command, " ");
	}
    
    
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
//            final ExecutionContext exec) throws Exception {
//
//    	String outfile = m_infile.getStringValue().replace(".bam","."+m_filesuffix.getStringValue());
//    	
//    	String cmd = "java -jar -Xmx"+m_memory.getIntValue()+"G " + m_path2gatk.getStringValue();
//    	cmd += " -T DepthOfCoverage";
//    	cmd += " -R "+m_refgenome.getStringValue();
//    	cmd += " -o "+outfile;
//    	cmd += " -I "+m_infile.getStringValue();
//    	
//    	if(!m_path2bed.getStringValue().equals("")){
//    		cmd += " -L "+m_path2bed.getStringValue();
//    	}
//    	
//    	cmd += " -nt "+m_threads.getIntValue();
//    	
//    	if(!m_extrafilters.getStringValue().equals("")){
//        	String[] filters = m_extrafilters.getStringValue().split(",");
//        	for(String filter : filters){
//        		cmd += " --read_filter "+filter;
//        	}
//    	}
//
//    	DepthOfCoverageNodeModel.logger.info("Running GATK DepthOfCoverage...");
//    	DepthOfCoverageNodeModel.logger.info("Log files can be found in "+outfile+".out.log and "+outfile+".err.log");
//    	
//		Executor.executeCommand(new String[]{cmd}, exec, null, DepthOfCoverageNodeModel.logger, outfile+".out.log", outfile+".err.log", null);
//		System.out.println("GATK finished, processing outfiles");
//    	DepthOfCoverage.processCoverageFile(outfile);
//		
//        return null;
//    }
    

	@Override
	protected String getCommandWalker() {
		return "DepthOfCoverage";
	}

	@Override
	protected File getLockFile() {
		return new File(this.LOCKFILE);
	}

	@Override
	protected String getOutfile() {
		return this.OUTFILE;
	}

	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_extrafilters.saveSettingsTo(settings);
        m_filesuffix.saveSettingsTo(settings);
        m_infile.saveSettingsTo(settings);
        m_path2bed.saveSettingsTo(settings);
        m_bed_file_checkbox.saveSettingsTo(settings);
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_extrafilters.loadSettingsFrom(settings);
        m_filesuffix.loadSettingsFrom(settings);
        m_infile.loadSettingsFrom(settings);
        m_path2bed.loadSettingsFrom(settings);
        m_bed_file_checkbox.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_extrafilters.validateSettings(settings);
        m_filesuffix.validateSettings(settings);
        m_infile.validateSettings(settings);
        m_path2bed.validateSettings(settings);
        m_bed_file_checkbox.validateSettings(settings);
	}
}

