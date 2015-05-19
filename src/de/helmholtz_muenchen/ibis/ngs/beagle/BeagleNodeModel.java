package de.helmholtz_muenchen.ibis.ngs.beagle;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission.GATKPhaseByTransmissionNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;


/**
 * This is the model implementation of Beagle.
 * 
 *
 * @author Tanzeem
 */
public class BeagleNodeModel extends NodeModel {
    
	static boolean optionalPort = false;

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_BEAGLE_PATH = "BEAGLE_PATH";
	public static final String CFGKEY_REF_VCF = "REFVCF";
	public static final String CFGKEY_INFILE = "INFILE";
	public static final String CFGKEY_PED_FILE = "PED";
    static final String CFGKEY_PARAM="parameter";
    static final String[] GENOTYPE_PARAMS={"Genotype", "Genotype probabilities", "Both"};
    private final SettingsModelString m_param = new SettingsModelString(CFGKEY_PARAM, "");
    private static final String INPUT_COL_NAME = "Path2VCFFile";
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_BEAGLE = new SettingsModelString(CFGKEY_BEAGLE_PATH, "");
    private final SettingsModelString m_INFILE = new SettingsModelString(CFGKEY_INFILE, "");
    private final SettingsModelString m_REF_VCF = new SettingsModelString(CFGKEY_REF_VCF, "");
    private final SettingsModelString m_PED_FILE = new SettingsModelString(CFGKEY_PED_FILE, "");
    
//	private static final String REFERENCE = "Reference";

	//The Output Col Names
	public static final String OUT_COL1 = "PHASED VARIANTS";
		
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKPhaseByTransmissionNodeModel.class); //not used yet
		
	/**
	 * Command line syntax: java -jar beagle.jar [arguments]

data input/output parameters ...
  gt=<VCF file: use GT field>                        (optional)
  gl=<VCF file: use GL/PL field>                     (optional)
  gtgl=<VCF file: use GT and GL/PL fields>           (optional)
  ref=<VCF file with phased genotypes>               (optional)
  out=<output file prefix>                           (required)
  excludesamples=<file with 1 sample ID per line>    (optional)
  excludemarkers=<file with 1 marker ID per line>    (optional)
  ped=<linkage format pedigree file>                 (optional)
  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)
  maxlr=<max GL/PL likelihood ratio>                 (default=5000)

algorithm parameters ...
  nthreads=<number of threads>                       (default=1)
  window=<markers per window>                        (default=50000)
  overlap=<overlap between windows>                  (default=3000)
  gprobs=<print GP field (true/false)>               (default=true)
  impute=<impute ungenotyped variants (true/false)>  (default=true)
  usephase=<use phase in "gt" or "gtgl" file>        (default=false)
  singlescale=<model scale for singles>              (default=1.0)
  duoscale=<model scale for duos>                    (default=1.0)
  trioscale=<model scale trios>                      (default=1.0)
  burnin-its=<number of iterations>                  (default=5)
  phase-its=<number of iterations>                   (default=5)
  impute-its=<number of iterations>                  (default=5)
  seed=<random seed>                                 (default=-99999)

IBD parameters ...
  ibd=<perform IBD detection (true/false)>           (default=false)
  ibdlod=<min LOD score for reporting IBD>           (default=3.0)
  ibdscale=<model scale factor for Refined IBD>      (default: data-dependent)
  ibdtrim=<markers at each segment end>              (default=40)

##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">
##FORMAT=<ID=PL,Number=G,Type=Float,Description="Phred-scaled Genotype Likelihoods">
	 */

    /**
     * Constructor for the node model.
     */
    protected BeagleNodeModel() {
    
        // TODO one incoming port and one outgoing port is assumed
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String inFile = "", outFile = "", refvcf = "";
    	
    	if (optionalPort) {
    		inFile = inData[0].iterator().next().getCell(0).toString(); //path2VCFFile     		
    	}
    	else {
        	inFile = m_INFILE.getStringValue();
    	}
    	
    	//checks if all files are still available
        if(!Files.exists(Paths.get(inFile))){ 
        	throw new Exception("file: "+ inFile+" does not exist");
        }
                
        LOGGER.info("Using file "+inFile);
        
        outFile = inFile.replaceAll(".vcf", "_beagle_phased");
    	refvcf = m_REF_VCF.getStringValue();
        String param=m_param.getStringValue();
        System.out.println("PARAM "+param);
      //check if tool is specified and node has been properly configured
        if(param.equals("")){
        	throw new Exception("You have to select the genotype parameter before running it!!!");
        }
//        LOGGER.info("Starting with the parame "+param+"...");

        String g = "";
        if (param.equals(BeagleNodeModel.GENOTYPE_PARAMS[0]))
        	g = "gt";
        if (param.equals(BeagleNodeModel.GENOTYPE_PARAMS[1]))
        	g = "gl";
        if (param.equals(BeagleNodeModel.GENOTYPE_PARAMS[2]))
        	g = "gtgl";

		
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-jar "+m_BEAGLE.getStringValue());
    	command.add(g + "=" + inFile);
    	command.add("out=" + outFile);
    	command.add("ped="+m_PED_FILE.getStringValue());

    	if(!refvcf.equals(""))
    		command.add("ref="+refvcf);
    	

    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
    	
    	/**
    	 * unzipping the beagle output 
    	 */
    	command = new ArrayList<String>(2);
    	if (Files.exists(Paths.get(outFile+".vcf")))
    		Files.delete(Paths.get(outFile+".vcf"));
    		
    	command.add("sh /home/ibis/tanzeem.haque/Documents/Scripts/Sequenciator/beagleUnzip.sh");
    	command.add(outFile+".vcf.gz");
        	
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);

    	
//    	java -jar $pathGATK/GenomeAnalysisTK.jar -R $absPath/$ref -T PhaseByTransmission 
//    	-V $absPath/$input -ped $absPath/family.ped -o $absPath/$outputDir/$outputVCF 
//    	-PF $absPath/$outputDir/$outputRT
    	
    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	outFile = outFile+".vcf";

    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outFile)};
    	
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
        // TODO Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
    	/**
    	 * Check OptionalInputPort
    	 */
		try{
//			inSpecs[0].getColumnNames();
			String[] colNames = inSpecs[0].getColumnNames();
			optionalPort=true;
			
//			for (String s : colNames)

			if (colNames.length > 1 && colNames[0].equals(INPUT_COL_NAME)) {
				m_INFILE.setEnabled(false);
			}
			else
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			

    	}catch(NullPointerException npe){
    		m_INFILE.setEnabled(true);
			optionalPort=false;
    	
    	}
		
    	if(m_REF_VCF.getStringValue().equals(""))
    		m_REF_VCF.setEnabled(false);
    	else
    		m_REF_VCF.setEnabled(true);

		DataColumnSpec[] allColSpecs = new DataColumnSpec[1];
    	allColSpecs[0] = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();

    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
//		System.out.println(optionalPort);
        return new DataTableSpec[]{outputSpec};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        // TODO save user settings to the config object.
    	m_param.saveSettingsTo(settings);
		m_BEAGLE.saveSettingsTo(settings);
		m_INFILE.saveSettingsTo(settings);
		m_PED_FILE.saveSettingsTo(settings);
		m_REF_VCF.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	m_param.loadSettingsFrom(settings);
    	m_BEAGLE.loadSettingsFrom(settings);
		m_INFILE.loadSettingsFrom(settings);
		m_PED_FILE.loadSettingsFrom(settings);
		m_REF_VCF.loadSettingsFrom(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
         
    	m_param.validateSettings(settings);
    	m_BEAGLE.validateSettings(settings);
		m_INFILE.validateSettings(settings);
		m_PED_FILE.validateSettings(settings);
		m_REF_VCF.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // TODO load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // TODO save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }

}

