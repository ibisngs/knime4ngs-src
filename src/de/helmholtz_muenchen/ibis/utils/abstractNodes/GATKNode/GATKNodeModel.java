package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortType;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

public abstract class GATKNodeModel extends HTExecutorNodeModel{

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_GATK_MEM = "GATK_MEM";
	public static final String CFGKEY_PATH2BED = "bed";
	public static final String CFGKEY_BED_FILE_CHECKBOX = "BED_FILE_CHECKBOX";
	public static final String CFGKEY_OPT_FLAGS ="OPT_FLAGS";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_path2bed = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_PATH2BED,"",true);
    private final SettingsModelBoolean m_bed_file_checkbox = new SettingsModelBoolean(GATKNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
    private final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(CFGKEY_OPT_FLAGS,"",false);

	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "OUTFILE";
	
	private boolean outtable = true;
    
    /**
     * Constructor for the node model.
     */
    protected GATKNodeModel(PortType[] INPORTS, PortType[] OUTPORTS) {
        super(INPORTS, OUTPORTS);
        if(OUTPORTS.length==0) {
        	outtable = false;
        }
    }
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-Xmx"+m_GATK_MEM.getIntValue()+"G");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T "+getCommandWalker());
    	command.add("-R "+m_REF_GENOME.getStringValue());    	
    	command.add(getCommandParameters(inData));
   
    	String OUTFILE = getOutfile(); 	
    	command.add("-o "+OUTFILE);
    	
    	if(m_bed_file_checkbox.getBooleanValue()){
			if(m_path2bed.getStringValue().equals("") || Files.notExists(Paths.get(m_path2bed.getStringValue()))) {
				setWarningMessage("Specify valid bed file!");
			} else {
				command.add("-L "+m_path2bed.getStringValue());
			}
    	}
    	
    	command.add(" "+m_OPT_FLAGS.getStringValue());
    	
    	LOGGER.info(StringUtils.join(command, " "));
     	super.executeCommand(new String[]{StringUtils.join(command, " ")},exec, getLockFile(),OUTFILE+".stdOut",OUTFILE+".stdErr");
     	
     	
    	/**
    	 * OUTPUT
    	 */
     	if(!outtable) {
     		return null;
     	}
     	
     	//Table1
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1_TABLE1, getOutColType()).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable1 = cont.getTable();


        return new BufferedDataTable[]{outTable1};
 
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
    	
    	if(!checkInputCellType(inSpecs)) {
    		new InvalidSettingsException("This node seems to be incompatible with the precedent node!");
    	}
    	
    	String gatk_warning = CheckUtils.checkSourceFile(m_GATK.getStringValue());
    	if(gatk_warning != null) {
    		setWarningMessage(gatk_warning);
    	}
    	
    	String ref_warning = CheckUtils.checkSourceFile(m_REF_GENOME.getStringValue());
    	if(ref_warning != null) {
    		setWarningMessage(ref_warning);
    	}
    	
		if (!outtable) {
			return null;
		}
		
	

		DataTableSpec outSpecTable1 = new DataTableSpec(
				new DataColumnSpec[] { new DataColumnSpecCreator(
						OUT_COL1_TABLE1, getOutColType()).createSpec() });
		return new DataTableSpec[] { outSpecTable1 };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
   	 	m_GATK.saveSettingsTo(settings);
   	 	m_REF_GENOME.saveSettingsTo(settings);
   	 	m_GATK_MEM.saveSettingsTo(settings);
   	 	m_path2bed.saveSettingsTo(settings);
   	 	m_bed_file_checkbox.saveSettingsTo(settings);
   	 	m_OPT_FLAGS.saveSettingsTo(settings);
   	 	saveExtraSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
		m_GATK.loadSettingsFrom(settings);
		m_REF_GENOME.loadSettingsFrom(settings);
		m_GATK_MEM.loadSettingsFrom(settings);
		m_path2bed.loadSettingsFrom(settings);
		m_bed_file_checkbox.loadSettingsFrom(settings);
		m_OPT_FLAGS.loadSettingsFrom(settings);
		loadExtraValidatedSettingsFrom(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
		m_GATK.validateSettings(settings);
		m_REF_GENOME.validateSettings(settings);
		m_GATK_MEM.validateSettings(settings);
		m_path2bed.validateSettings(settings);
		m_bed_file_checkbox.validateSettings(settings);
		m_OPT_FLAGS.validateSettings(settings);
		validateExtraSettings(settings);
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

    /****************************** ABSTRACT METHODS **********************************/
    /**
     * Provides the node specific filter settings
     * @return
     */
    protected abstract String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException;
    	
    /**
     * Provides the GATK Walker
     * @return
     */
    protected abstract String getCommandWalker();
    protected abstract File getLockFile();
    protected abstract String getOutfile();
    protected abstract boolean checkInputCellType(DataTableSpec[] inSpecs);
    protected abstract DataType getOutColType();
    
    protected abstract void saveExtraSettingsTo(final NodeSettingsWO settings);
    protected abstract void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException;
    protected abstract void validateExtraSettings(final NodeSettingsRO settings) throws InvalidSettingsException;
    
}
