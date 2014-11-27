package de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants;

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

/**
 * This is the model implementation of GATKSelectVariants.
 * 
 *
 * @author Maximilian Hastreiter
 */
public abstract class SelectVariantsNodeModel extends NodeModel {
    
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_VCFCOLUMN = "VCFCOLUMN";
	public static final String CFGKEY_FILTERSTRING = "FILTERSTRING";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_VCFCOLUMN = new SettingsModelIntegerBounded(CFGKEY_VCFCOLUMN, 1, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_FILTERSTRING = new SettingsModelString(CFGKEY_FILTERSTRING, "");

	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(SelectVariantsNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "FilteredVCF";
    
    /**
     * Constructor for the node model.
     */
    protected SelectVariantsNodeModel(int INPORTS, int OUTPORTS) {
        super(INPORTS, OUTPORTS);
    }
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	/**
    	 * Get Input Data
    	 */
    	String INVCF = inData[0].iterator().next().getCell(m_VCFCOLUMN.getIntValue()-1).toString();
    	if(!INVCF.contains(".vcf")){
    		throw new InvalidSettingsException("Infile "+INVCF+" seems to be in the wrong format. VCF File required!");
    	}
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T SelectVariants");
    	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-V "+INVCF);
    	
    	command.add(getCommandParameters());
    	
    	String OUTFILE = INVCF.replace(".vcf", getOutfileSuffix());
    	
    	command.add("-o "+OUTFILE);
    	
    	System.out.println(StringUtils.join(command, " "));
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
     	
    	
     	
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
    		
    		if(m_VCFCOLUMN.getIntValue()>inSpecs[0].getNumColumns()){
    			throw new InvalidSettingsException("Selected column "+m_VCFCOLUMN.getIntValue()+" is not available. inData has only "+inSpecs[0].getNumColumns()+" columns!");
    		}
    		
    		
    		
        return new DataTableSpec[]{outSpecTable1};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
   	 m_GATK.saveSettingsTo(settings);
   	 m_REF_GENOME.saveSettingsTo(settings);
   	 m_VCFCOLUMN.saveSettingsTo(settings);
   	 m_FILTERSTRING.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
      	 m_GATK.loadSettingsFrom(settings);
       	 m_REF_GENOME.loadSettingsFrom(settings);
       	 m_VCFCOLUMN.loadSettingsFrom(settings);
       	 m_FILTERSTRING.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
      	 m_GATK.validateSettings(settings);
       	 m_REF_GENOME.validateSettings(settings);
       	 m_VCFCOLUMN.validateSettings(settings);
       	 m_FILTERSTRING.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

    protected SettingsModelString getFILTERSTRINGModel(){
    	return m_FILTERSTRING;
    }
    
    /****************************** ABSTRACT METHODS **********************************/
    /**
     * Provides the node specific filter settings
     * @return
     */
    protected abstract String getCommandParameters();
    
    /**
     * Provides the outfile suffix
     * @return
     */
    protected abstract String getOutfileSuffix();
    
}

