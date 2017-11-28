package de.helmholtz_muenchen.ibis.ngs.smartPhase;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.trimgalore.TrimGaloreNodeModel;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

/**
 * This is the model implementation of SmartPhase. KNIME integration of the
 * smartPhase program that phases filtered rare variants in desired genomic
 * regions to determine possible compound heterozygosity.
 *
 * @author Paul Hager
 */
public class SmartPhaseNodeModel extends HTExecutorNodeModel {
	
	// keys for SettingsModels
    protected static final String CFGKEY_FILE_LIST 			= "FileList";
    protected static final String CFGKEY_FILE_DIR 			= "FileDir";
    protected static final String CFGKEY_FILE_FILE 			= "FileFile";
    protected static final String CFGKEY_FILE_LIST_DISPLAY 	= "FileListDisplay";
    protected static final String CFGKEY_REGEX			 	= "FilenameRegex";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_FILE_DIR = "-";
    protected static final String DEFAULT_FILE_FILE = "-";
    protected static final String DEFAULT_REGEX = ".*";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FILE_DIR 					= new SettingsModelString(CFGKEY_FILE_DIR, DEFAULT_FILE_DIR);
    private final SettingsModelString SET_FILE_FILE 				= new SettingsModelString(CFGKEY_FILE_FILE, DEFAULT_FILE_FILE);
    private final SettingsModelString SET_NAME_REGEX				= new SettingsModelString(CFGKEY_REGEX, DEFAULT_REGEX);
    
    // storage 
    private final HashSet<String> FILES		= new HashSet<String>();
//    private boolean hasConfigureOpendOnce 	= false; // true, if configure was opend once
	
	// the logger instance
    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(TrimGaloreNodeModel.class);

    // Node specific params
    public static final String CFGKEY_SMARTPHASE = "Path2SmartPhase";
    public static final String CFGKEY_GENE_REGIONS = "Path2BEDFile";
    public static final String CFGKEY_FILT_VARS = "Path2FilteredVariantsFile";
    public static final String CFGKEY_ALL_VARS = "Path2AllVarsVCF";
    public static final String CFGKEY_READS = "Path2AllReadFiles";
    public static final String CFGKEY_MIN_MAPQ = "MinimumMappingQuality";
    public static final String CFGKEY_PATIENT_ID = "PatientID";
    public static final String CFGKEY_OUTPUT = "Path2Output";
    public static final String CFGKEY_TRIO = "TrioBoolean";
    public static final String CFGKEY_PED = "Path2Ped";
    
    private final SettingsModelString m_smartphase = 
			new SettingsModelString(CFGKEY_SMARTPHASE,"");
    
    private final SettingsModelString m_gene_regions = 
			new SettingsModelString(CFGKEY_GENE_REGIONS,"");
    
    private final SettingsModelString m_filt_vars = 
			new SettingsModelString(CFGKEY_FILT_VARS,"");
    
    private final SettingsModelString m_all_vars = 
			new SettingsModelString(CFGKEY_ALL_VARS,"");
    
    private final SettingsModelString m_patient_id = 
			new SettingsModelString(CFGKEY_PATIENT_ID,"");
    
    private final SettingsModelString m_output = 
			new SettingsModelString(CFGKEY_OUTPUT,"");
    
    private final SettingsModelString m_ped = 
			new SettingsModelString(CFGKEY_PED,"");
    
    private final SettingsModelBoolean m_trio =
    		new SettingsModelBoolean(CFGKEY_TRIO, false);
    
    private final SettingsModelString m_mapq =
    		new SettingsModelString(CFGKEY_MIN_MAPQ, "");
    
    //The Output Col Names
  	public static final String OUT_SMARTPHASE_PRED = "Path2OutputFilePredictions";
  	public static final String OUT_SMARTPHASE_STAT = "Path2OutputFileStatistics";
  	
  	public static final String OUTPUT_NAME_BAM_FILES 	= "BAMFiles";
	public static final String FILTERTYPE_NAME 			= "BAM";
	public static final String FILENAME_END_REGEX 		= "\\.bam$";
    
	/**
	 * Constructor for the node model.
	 */
	protected SmartPhaseNodeModel() {

		super(1, 1, 1);
		
		addSetting(m_smartphase);
		addSetting(m_gene_regions);
		addSetting(m_filt_vars);
		addSetting(m_all_vars);
		addSetting(m_patient_id);
		addSetting(m_output);
		addSetting(m_ped);
		addSetting(m_trio);
		addSetting(m_mapq);
		addSetting(SET_FILE_DIR);
		addSetting(SET_FILE_FILE);
		addSetting(SET_NAME_REGEX);
		
		addPrefPageSetting(m_smartphase, IBISKNIMENodesPlugin.SMARTPHASE);
		
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec)
			throws Exception {
		
		int colIndex1 = 0;
    	if(SmartPhaseNodeDialog.getUseMainInputColBool()){
    		colIndex1 = inData[0].getDataTableSpec().findColumnIndex(SmartPhaseNodeDialog.getMainInputCol1());
    	}
    	String allVarsFilePath = inData[0].iterator().next().getCell(colIndex1).toString();
		
		String smartPhasePath = m_smartphase.getStringValue();
		String genomicRegionsFilePath = m_gene_regions.getStringValue();
		String filtVarsFilePath = m_filt_vars.getStringValue();
		String patientID = m_patient_id.getStringValue();
		String outputPath = m_output.getStringValue();
		String pedFilPath = m_ped.getStringValue();
		boolean trio = m_trio.getBooleanValue();
		String minMapQs = m_mapq.getStringValue();
		
		// Beautify outfolder path
    	outputPath = outputPath.trim();
    	if(!outputPath.equals("") && !outputPath.endsWith(System.getProperty("file.separator"))){
    		outputPath += System.getProperty("file.separator");
    	}
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	String outFilePred = "";
    	String outFileStat = "";
    	
    	if(!outputPath.equals("")){
    		outFilePred = outputPath+patientID+".smartPase.predictions.out";
    		outFileStat = outputPath+patientID+".smartPase.statistics.out";
    	}
    	
    	File lockFile = new File(outFilePred.substring(0,outFilePred.lastIndexOf(".")) + ".sp" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	ArrayList<String> cmd = new ArrayList<String>();
    	cmd.add("java");
    	cmd.add("-jar");
    	cmd.add(IO.processFilePath(smartPhasePath));
    	
    	cmd.add("-g="+genomicRegionsFilePath);
    	
    	cmd.add("-f="+filtVarsFilePath);
    	
    	cmd.add("-a="+allVarsFilePath);
    	
    	cmd.add("-p="+patientID);
    	
    	cmd.add("-o="+outFilePred);
    	
    	String BAMFiles = "";
    	for(String fileName : FILES){
    		BAMFiles += fileName;
    		BAMFiles += ",";
    	}
    	if(!BAMFiles.equals("")){
    		BAMFiles.substring(0, BAMFiles.length()-1);
    		cmd.add("-r="+BAMFiles);
    		cmd.add("-m="+minMapQs);
    	}
    	
    	if(trio){
    		cmd.add("-t");	
    		cmd.add("-d="+pedFilPath);
    	}
    	
    	// Parse options into command string array - ensures spaces in path will be correctly interpreted
    	String[] command = new String[cmd.size()];
    	for(int indx=0; indx < cmd.size(); indx++){
    		command[indx] = cmd.get(indx);
    	}
    	
    	LOGGER.info("CMD: "+command.toString());
    	
    	LOGGER.info("Running smartPhase...");
		LOGGER.info("Log files can be found in "+outFilePred+".stdOut and "+outFilePred+".stdErr");
		//command = new String[0];
		//command[0] = "ll";
		for(String s : command){
			System.out.println(s);
		}
		super.executeCommand(command, outFilePred, exec, new String[]{}, lockFile, outFilePred+".stdOut", outFilePred+".stdErr", null, null, null);
    	
		
		BufferedDataContainer cont;
	    FileCell[] c;
	    	
	   
	    cont = exec.createDataContainer(createSpecs());
	        
	    
	    c = new FileCell[]{
	   			FileCellFactory.create(outFilePred),
	   			FileCellFactory.create(outFileStat)};
	   
	    cont.addRowToTable(new DefaultRow("Row0", c));
	  
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	

		// TODO: Return a BufferedDataTable for each output port
		return new BufferedDataTable[] {outTable};
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

		CompatibilityChecker CC = new CompatibilityChecker();
    	if(CC.getWarningStatus()){
    		setWarningMessage(CC.getWarningMessages());
    	}
    	
    	String smartPhase = m_smartphase.getStringValue();
    	if(CompatibilityChecker.inputFileNotOk(smartPhase) || !smartPhase.endsWith(".jar")) {
    		throw new InvalidSettingsException("Provided smartPhase.jar does not exist or is invalid!");
    	}
    	
    	String geneRegions = m_gene_regions.getStringValue();
    	if(CompatibilityChecker.inputFileNotOk(geneRegions) || !geneRegions.endsWith(".bed")) {
    		throw new InvalidSettingsException("Provided genomic regions bed file does not exist or is invalid! Please provide a valid .bed file.");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_filt_vars.getStringValue())) {
    		throw new InvalidSettingsException("Provided filtered variants file does not exist or is invalid!");
    	}
    	
    	String allVars = m_all_vars.getStringValue();
    	if(CompatibilityChecker.inputFileNotOk(allVars) || (!allVars.endsWith(".vcf") && !allVars.endsWith(".vcf.gz"))) {
    		throw new InvalidSettingsException("Provided all variants file does not exist or is invalid! Please provide a valid .vcf file.");
    	}
    	
    	if(m_trio.getBooleanValue()){
    		String ped = m_ped.getStringValue();
        	if(CompatibilityChecker.inputFileNotOk(ped) || !ped.endsWith(".ped")) {
        		throw new InvalidSettingsException("Provided pedigree file does not exist or is invalid! Please provide a valid .ped file.");
        	}
    	}
    	
		
		return new DataTableSpec[] { createSpecs() };
	}
	
	/**
     * Create Tablespecs
     * @return
     */
    private DataTableSpec createSpecs(){
    	DataTableSpec out;
    	
    	out = new DataTableSpec(
            	new DataColumnSpec[]{
            			new DataColumnSpecCreator(OUT_SMARTPHASE_PRED, FastQCell.TYPE).createSpec(),
            			new DataColumnSpecCreator(OUT_SMARTPHASE_STAT, FastQCell.TYPE).createSpec()});
    	    	  	
    	return out;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	this.FILES.clear();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	// clean the old data
    	FILES.clear();
    	// check, if data is set
        if (settings.containsKey(SmartPhaseNodeModel.CFGKEY_FILE_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(SmartPhaseNodeModel.CFGKEY_FILE_LIST))
					this.FILES.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }
    
    @Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
 
    	settings.addStringArray(SmartPhaseNodeModel.CFGKEY_FILE_LIST, FILES.toArray(new String[FILES.size()]));
	}

	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		super.validateSettings(settings);
	}
	
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {	
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
	}

}
