package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

/**
 * This is the model implementation of GATKPhaseByTransmission.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKPhaseByTransmissionNodeModel extends NodeModel {
    
//	my $com = "java -Xmx8g -jar ".$TOOL_PATH.$GATK_VERSION."/GenomeAnalysisTK.jar -T PhaseByTransmission
//	-R $REF_GENOME -V $INFILE -o $outfile -ped $PED";
	
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_INFILE = "INFILE";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_PED_FILE = "PED";
	public static final String CFGKEY_JAVAMEMORY = "gatkmemory";
	public static final String CFGKEY_DENOVOPRIOR = "deNovoPrior";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "---");
    private final SettingsModelString m_INFILE = new SettingsModelString(CFGKEY_INFILE, "---");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "---");
    private final SettingsModelString m_PED_FILE = new SettingsModelString(CFGKEY_PED_FILE, "---");

    
    public static final int DEF_NUM_JAVAMEMORY=8;
    public static final int MIN_NUM_JAVAMEMORY=1;
    public static final int MAX_NUM_JAVAMEMORY=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_GATK_JAVA_MEMORY = new SettingsModelIntegerBounded(CFGKEY_JAVAMEMORY, DEF_NUM_JAVAMEMORY, MIN_NUM_JAVAMEMORY, MAX_NUM_JAVAMEMORY);
    
    private final SettingsModelString m_DENOVO_PRIOR = new SettingsModelString(CFGKEY_DENOVOPRIOR, "1.0E-8");
    
	//The Output Col Names
	public static final String OUT_COL1 = "PHASED VARIANTS";
	
    /**
     * Constructor for the node model.
     */
    protected GATKPhaseByTransmissionNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-Xmx"+m_GATK_JAVA_MEMORY.getIntValue()+"g -jar "+m_GATK.getStringValue());
    	command.add("-T PhaseByTransmission");
    	
    	command.add("--DeNovoPrior "+m_DENOVO_PRIOR.getStringValue());
    	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-V "+m_INFILE.getStringValue());
    	
    	String OUTFILE = m_INFILE.getStringValue().replaceAll(".vcf", "_phased.vcf");
    	command.add("-o "+OUTFILE);
    	
    	command.add("-ped "+m_PED_FILE.getStringValue());
    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();

        return new BufferedDataTable[]{outTable};
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

        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
		   	 m_GATK.saveSettingsTo(settings);
		   	 m_INFILE.saveSettingsTo(settings);
		   	 m_REF_GENOME.saveSettingsTo(settings);
		   	 m_PED_FILE.saveSettingsTo(settings);
		   	 m_GATK_JAVA_MEMORY.saveSettingsTo(settings);
		   	 m_DENOVO_PRIOR.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_GATK.loadSettingsFrom(settings);
		   	 m_INFILE.loadSettingsFrom(settings);
		   	 m_REF_GENOME.loadSettingsFrom(settings);
		   	 m_PED_FILE.loadSettingsFrom(settings);
		   	 m_GATK_JAVA_MEMORY.loadSettingsFrom(settings);
		   	 m_DENOVO_PRIOR.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_GATK.validateSettings(settings);
		   	 m_INFILE.validateSettings(settings);
		   	 m_REF_GENOME.validateSettings(settings);
		   	 m_PED_FILE.validateSettings(settings);
		     m_GATK_JAVA_MEMORY.validateSettings(settings);
		     m_DENOVO_PRIOR.validateSettings(settings);
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

}

