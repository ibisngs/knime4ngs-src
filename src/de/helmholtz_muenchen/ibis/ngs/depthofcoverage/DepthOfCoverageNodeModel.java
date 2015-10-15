package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
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
	public static final String CFGKEY_INFILE = "infile";
	public static final String CFGKEY_EXTRAFILTERS = "extrafilters";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	/**
	 * Node Models
	 */
	private final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
	private final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"DoC");
	
	private String OUTFILE, LOCKFILE;
	private int bam_index;

    /**
     * Constructor for the node model.
     */
    protected DepthOfCoverageNodeModel() {
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

    @Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {

    	/**
    	 * Check INFILE
    	 */
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(bam_index).toString();
    		if(!INFILE.endsWith(".bam")){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to BAM infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to BAM infile but it is empty.");
    	}
		
    	ArrayList<String> command = new ArrayList<String>();
		
		String fileSuffix = m_filesuffix.getStringValue();
		if(fileSuffix.equals("")) {
			fileSuffix = "DoC";
		}
		
		this.OUTFILE = IO.replaceFileExtension(INFILE, fileSuffix);
		this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
		command.add("-I "+INFILE);

		
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
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_extrafilters.loadSettingsFrom(settings);
        m_filesuffix.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_extrafilters.validateSettings(settings);
        m_filesuffix.validateSettings(settings);
	}

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		bam_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("BAMCell")) {
    			bam_index = i;
    		}
    	}
		return (bam_index>-1);
	}

	@Override
	protected DataType getOutColType() {
		return FileCell.TYPE;
	}
}

