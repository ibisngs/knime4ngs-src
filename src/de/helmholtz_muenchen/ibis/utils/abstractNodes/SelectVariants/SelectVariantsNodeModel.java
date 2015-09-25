package de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKSelectVariants.
 * 
 *
 * @author Maximilian Hastreiter
 */
public abstract class SelectVariantsNodeModel extends GATKNodeModel {
    
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_VCFCOLUMN = "VCFCOLUMN";
	public static final String CFGKEY_FILTERSTRING = "FILTERSTRING";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelIntegerBounded m_VCFCOLUMN = new SettingsModelIntegerBounded(CFGKEY_VCFCOLUMN, 1, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_FILTERSTRING = new SettingsModelString(CFGKEY_FILTERSTRING, "");
    
	private String OUTFILE, LOCKFILE;
	
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "FilteredVCF";
    
    /**
     * Constructor for the node model.
     */
    protected SelectVariantsNodeModel(int INPORTS, int OUTPORTS) {
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

    protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {
    	ArrayList<String> command = new ArrayList<String>();
    	
    	String INVCF = inData[0].iterator().next().getCell(m_VCFCOLUMN.getIntValue()-1).toString();
    	if(!INVCF.contains(".vcf")){
    		throw new InvalidSettingsException("Infile "+INVCF+" seems to be in the wrong format. VCF File required!");
    	}
    	command.add("-V "+INVCF);
    	command.add(getCommandParameters());
    	this.OUTFILE = IO.replaceFileExtension(INVCF, getOutfileSuffix());
		this.LOCKFILE = IO.replaceFileExtension(INVCF, SuccessfulRunChecker.LOCK_ENDING);
    	
    	return StringUtils.join(command, " ");
    }
    
    protected File getLockFile() {
    	return new File(this.LOCKFILE);
    }
    
    
    protected String getOutfile() {
    	return this.OUTFILE;
    }

//    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
//            final ExecutionContext exec) throws Exception {
//
//    	/**
//    	 * Get Input Data
//    	 */
//    	String INVCF = inData[0].iterator().next().getCell(m_VCFCOLUMN.getIntValue()-1).toString();
//    	if(!INVCF.contains(".vcf")){
//    		throw new InvalidSettingsException("Infile "+INVCF+" seems to be in the wrong format. VCF File required!");
//    	}
//    	
//    	ArrayList<String> command = new ArrayList<String>();
//    	
//    	command.add("java");
//    	command.add("-jar "+m_GATK.getStringValue());
//    	command.add("-T SelectVariants");
//    	
//    	command.add("-R "+m_REF_GENOME.getStringValue());
//    	command.add("-V "+INVCF);
//    	
//    	command.add(getCommandParameters());
//    	
//    	String OUTFILE = INVCF.replace(".vcf", getOutfileSuffix());
//    	
//    	command.add("-o "+OUTFILE);
//    	
//    	System.out.println(StringUtils.join(command, " "));
//     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
//     	
//    	
//     	
//    	/**
//    	 * OUTPUT
//    	 */
//     	//Table1
//    	BufferedDataContainer cont = exec.createDataContainer(
//    			new DataTableSpec(
//    			new DataColumnSpec[]{
//    					new DataColumnSpecCreator(OUT_COL1_TABLE1, FileCell.TYPE).createSpec()}));
//    	
//    	FileCell[] c = new FileCell[]{
//    			(FileCell) FileCellFactory.create(OUTFILE)};
//    	
//    	cont.addRowToTable(new DefaultRow("Row0",c));
//    	cont.close();
//    	BufferedDataTable outTable1 = cont.getTable();
//
//
//        return new BufferedDataTable[]{outTable1};
// 
//    }

    
//    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
//            throws InvalidSettingsException {
//    	
//    		DataTableSpec outSpecTable1 = new DataTableSpec(
//    										new DataColumnSpec[]{
//    										new DataColumnSpecCreator(OUT_COL1_TABLE1, FileCell.TYPE).createSpec()});
//    		
//    		if(m_VCFCOLUMN.getIntValue()>inSpecs[0].getNumColumns()){
//    			throw new InvalidSettingsException("Selected column "+m_VCFCOLUMN.getIntValue()+" is not available. inData has only "+inSpecs[0].getNumColumns()+" columns!");
//    		}
//    		
//    		
//    		
//        return new DataTableSpec[]{outSpecTable1};
//    }

    protected String getCommandWalker() {
    	return "SelectVariants";
    }

    protected void saveExtraSettingsTo(final NodeSettingsWO settings) {
   	 m_VCFCOLUMN.saveSettingsTo(settings);
   	 m_FILTERSTRING.saveSettingsTo(settings);
    }

    protected void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
       	 m_VCFCOLUMN.loadSettingsFrom(settings);
       	 m_FILTERSTRING.loadSettingsFrom(settings);
    }

    protected void validateExtraSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
       	 m_VCFCOLUMN.validateSettings(settings);
       	 m_FILTERSTRING.validateSettings(settings);
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

