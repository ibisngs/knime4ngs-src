package de.helmholtz_muenchen.ibis.ngs.gatkselectvariants;

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
public class GATKSelectVariantsNodeModel extends NodeModel {
    
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_VCFCOLUMN = "VCFCOLUMN";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_VCFCOLUMN = new SettingsModelIntegerBounded(CFGKEY_VCFCOLUMN, 1, 1, Integer.MAX_VALUE);
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKSelectVariantsNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "VCF_SNPs";
	public static final String OUT_COL1_TABLE2 = "VCF_Indels";
    
    /**
     * Constructor for the node model.
     */
    protected GATKSelectVariantsNodeModel() {
        super(1, 2);
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
    	
    	/**
    	 * Prepare Command for Indels
    	 */
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T SelectVariants");
    	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-V "+INVCF);
    	
    	String OUTFILE_INDELS = INVCF.replaceAll(".vcf", "_INDELS.vcf");
    	command.add("-o "+OUTFILE_INDELS);
    	command.add("-selectType INDEL");
    	
    	System.out.println(StringUtils.join(command, " "));
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
     	
    	/**
    	 * Prepare Command for SNPs
    	 */
    	command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T SelectVariants");
    	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-V "+INVCF);
    	
    	String OUTFILE_SNPS = INVCF.replaceAll(".vcf", "_SNPS.vcf");
    	command.add("-o "+OUTFILE_SNPS);
    	command.add("-selectType SNP");
    	
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
    			(FileCell) FileCellFactory.create(OUTFILE_SNPS)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable1 = cont.getTable();
    	
    	//Table2
    	cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1_TABLE2, FileCell.TYPE).createSpec()}));
    	
    	c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE_INDELS)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable2 = cont.getTable();

        return new BufferedDataTable[]{outTable1,outTable2};
 
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
    		DataTableSpec outSpecTable2 = new DataTableSpec(
					new DataColumnSpec[]{
					new DataColumnSpecCreator(OUT_COL1_TABLE2, FileCell.TYPE).createSpec()});

    		
    		if(m_VCFCOLUMN.getIntValue()>inSpecs[0].getNumColumns()){
    			throw new InvalidSettingsException("Selected column "+m_VCFCOLUMN.getIntValue()+" is not available. inData has only "+inSpecs[0].getNumColumns()+" columns!");
    		}
    		
    		
    		
        return new DataTableSpec[]{outSpecTable1,outSpecTable2};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
   	 m_GATK.saveSettingsTo(settings);
   	 m_REF_GENOME.saveSettingsTo(settings);
   	 m_VCFCOLUMN.saveSettingsTo(settings);
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

