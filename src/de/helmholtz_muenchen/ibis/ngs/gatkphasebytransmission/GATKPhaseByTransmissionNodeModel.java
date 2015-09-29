package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKPhaseByTransmission.
 * 
 *
 * @author Maximilian Hastreiter
 * modified by Tanzeem Haque
 * modified by Tim Jeske
 */
public class GATKPhaseByTransmissionNodeModel extends GATKNodeModel {

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_PED_FILE = "PED";
	public static final String CFGKEY_DENOVOPRIOR = "deNovoPrior";

	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_PED_FILE = new SettingsModelString(CFGKEY_PED_FILE, ""); 
    private final SettingsModelString m_DENOVO_PRIOR = new SettingsModelString(CFGKEY_DENOVOPRIOR, "1.0E-8");
    
	//The Output Col Names
	public static final String OUT_COL1 = "PhasedVCF";
		
	private String OUTFILE, LOCKFILE;
	
    /**
     * Constructor for the node model.
     */
    protected GATKPhaseByTransmissionNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData)
			throws InvalidSettingsException {
		/**
    	 * Check INFILE
    	 */
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(0).toString();
    		if(!INFILE.endsWith(".vcf")){
    			throw new InvalidSettingsException("First cell of input table has to be the path to VCF infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("First cell of input table has to be the path to VCF infile but it is empty.");
    	}
    	
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("-V "+INFILE);
    	
    	OUTFILE = IO.replaceFileExtension(INFILE, ".PhasedByTransmission.vcf");
    	LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
    	
    	command.add("-prior "+m_DENOVO_PRIOR.getStringValue());
    	command.add("-ped "+m_PED_FILE.getStringValue());
    	
    	return StringUtils.join(command, " ");
	}
	
	
    /**
     * {@inheritDoc}
     */
//    @Override
//    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
//            final ExecutionContext exec) throws Exception {
//    	
//    	String refGenome = "", inFile = "";
//    	
//    	if (optionalPort) {
//    		inFile = inData[0].iterator().next().getCell(0).toString(); //path2VCFFile
//    		refGenome = getAvailableFlowVariables().get(REFERENCE).getStringValue();
//     		
//    	}
//    	else {
//        	inFile = m_INFILE.getStringValue();
//        	refGenome = m_REF_GENOME.getStringValue();
//    	}
//    	
////    	System.out.println("infile "+inFile);
////		System.out.println("refgenome: "+refGenome);
//		
//    	ArrayList<String> command = new ArrayList<String>();
//    	
//    	command.add("java");
//    	command.add("-Xmx"+m_GATK_JAVA_MEMORY.getIntValue()+"g -jar "+m_GATK.getStringValue());
//    	command.add("-T PhaseByTransmission");
//    	
//    	command.add("-prior "+m_DENOVO_PRIOR.getStringValue());
//    	
//    	command.add("-R "+refGenome);
//    	command.add("-V "+inFile);
//    	
//    	String OUTFILE = inFile.replaceAll(".vcf", "_pbt_phased.vcf");
////		System.out.println("outfile: "+OUTFILE);
//
//    	command.add("-o "+OUTFILE);
//    	
//    	command.add("-ped "+m_PED_FILE.getStringValue());
//
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
////    	java -jar $pathGATK/GenomeAnalysisTK.jar -R $absPath/$ref -T PhaseByTransmission 
////    	-V $absPath/$input -ped $absPath/family.ped -o $absPath/$outputDir/$outputVCF 
////    	-PF $absPath/$outputDir/$outputRT
//    	
//    	
//    	/**
//    	 * OUTPUT
//    	 */
//    	BufferedDataContainer cont = exec.createDataContainer(
//    			new DataTableSpec(
//    			new DataColumnSpec[]{
//    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
//    	
//    	FileCell[] c = new FileCell[]{
//    			(FileCell) FileCellFactory.create(OUTFILE)};
//    	
//    	cont.addRowToTable(new DefaultRow("Row0",c));
//    	cont.close();
//    	BufferedDataTable outTable = cont.getTable();
//
//        return new BufferedDataTable[]{outTable};
//    }

    /**
     * {@inheritDoc}
     */
//    @Override
//    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
//            throws InvalidSettingsException {
//
//    	/**
//    	 * Check OptionalInputPort
//    	 */
//		try{
////			inSpecs[0].getColumnNames();
//			String[] colNames = inSpecs[0].getColumnNames();
//			optionalPort=true;
//			
////			for (String s : colNames)
//
//			if (colNames.length > 1 && colNames[0].equals(INPUT_COL_NAME)) {
//				m_INFILE.setEnabled(false);
//				if (getAvailableInputFlowVariables().containsKey(REFERENCE))
//					m_REF_GENOME.setEnabled(false);
////				System.out.println(colNames[0]);
//			}
//			else
//				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
//			
//
//    	}catch(NullPointerException npe){
//    		m_INFILE.setEnabled(true);
//    		m_REF_GENOME.setEnabled(true);
//			optionalPort=false;
//    	
//    	}
//		
//    	DataColumnSpec[] allColSpecs = new DataColumnSpec[1];
//    	allColSpecs[0] = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();
//
//    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
////		System.out.println(optionalPort);
//        return new DataTableSpec[]{outputSpec};
//        // TODO: generated method stub
//    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveExtraSettingsTo(final NodeSettingsWO settings) {
		   	 m_PED_FILE.saveSettingsTo(settings);
		   	 m_DENOVO_PRIOR.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_PED_FILE.loadSettingsFrom(settings);
		   	 m_DENOVO_PRIOR.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateExtraSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_PED_FILE.validateSettings(settings);
		     m_DENOVO_PRIOR.validateSettings(settings);
    } 

	@Override
	protected String getCommandWalker() {
		return "PhaseByTransmission";
	}

	@Override
	protected File getLockFile() {
		return new File(LOCKFILE);
	}

	@Override
	protected String getOutfile() {
		return OUTFILE;
	}
}

