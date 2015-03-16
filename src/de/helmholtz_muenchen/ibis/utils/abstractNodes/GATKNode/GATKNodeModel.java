package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

public abstract class GATKNodeModel extends NodeModel{

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_GATK_MEM = "GATK_MEM";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);


	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "OUTFILE";
    
    /**
     * Constructor for the node model.
     */
    protected GATKNodeModel(int INPORTS, int OUTPORTS) {
        super(INPORTS, OUTPORTS);
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
    	
    	System.out.println(StringUtils.join(command, " "));
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,null,LOGGER,OUTFILE+".stdOut",OUTFILE+".stdErr");
     	
    	
     	
    	/**
    	 * OUTPUT
    	 */
     	//Table1
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1_TABLE1, FileCell.TYPE).createSpec()}));
    	
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
    	
    		DataTableSpec outSpecTable1 = new DataTableSpec(
    										new DataColumnSpec[]{
    										new DataColumnSpecCreator(OUT_COL1_TABLE1, FileCell.TYPE).createSpec()});    		
        return new DataTableSpec[]{outSpecTable1};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
   	 	m_GATK.saveSettingsTo(settings);
   	 	m_REF_GENOME.saveSettingsTo(settings);
   	 	m_GATK_MEM.saveSettingsTo(settings);
   	 	saveExtraSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
      	 m_GATK.loadSettingsFrom(settings);
       	 m_REF_GENOME.loadSettingsFrom(settings);
       	 m_GATK_MEM.loadSettingsFrom(settings);
       	 loadExtraValidatedSettingsFrom(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
      	 m_GATK.validateSettings(settings);
       	 m_REF_GENOME.validateSettings(settings);
       	 m_GATK_MEM.validateSettings(settings);
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
    protected abstract String getCommandParameters(final BufferedDataTable[] inData);
    	
    /**
     * Provides the GATK Walker
     * @return
     */
    protected abstract String getCommandWalker();
    
    protected abstract String getOutfile();
    
    protected abstract void saveExtraSettingsTo(final NodeSettingsWO settings);
    protected abstract void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException;
    protected abstract void validateExtraSettings(final NodeSettingsRO settings) throws InvalidSettingsException;
    
}
