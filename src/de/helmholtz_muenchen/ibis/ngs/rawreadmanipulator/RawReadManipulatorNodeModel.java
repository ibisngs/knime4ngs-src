package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulator;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


/**
 * This is the model implementation of RawReadManipulator.
 * 
 *
 * @author Sebastian Kopetzky 
 * @author Maximilian Hastreiter
 */
public class RawReadManipulatorNodeModel extends HTExecutorNodeModel {
    

	//The Output Col Names
	public static final String OUT_COL1 = "Path2ReadFile1";
	public static final String OUT_COL2 = "Path2ReadFile2";
	
	public static final String CFGKEY_FILTERFILEEXISTS = "filterfileexists";
	public static final String CFGKEY_USEOTHERFILTERFILE = "useownfilterfile";
	public static final String CFGKEY_OTHERFILTERSETTINGSFILE = "otherfiltersettingsfile";

	public static final String CFGKEY_IFBARCODEFILE = "ifbarcodefile";
	public static final String CFGKEY_BARCODEFILE = "barcodefile";
	
	public static final String CFGKEY_REMOVEADAPTERS = "removeadapters";
	public static final String CFGKEY_ADAPTERS = "adapters";

	public static final String CFGKEY_DOTRIMPOLYAT = "dotrimpolyat";
	public static final String CFGKEY_TRIMPOLYAT = "trimpolyat";
	
	public static final String CFGKEY_LENGTHCUTOFF = "lengthcutoff";
	public static final String CFGKEY_MINLENGTH = "minlength";

	public static final String CFGKEY_PRESERVE = "preserve"; //if flag present, keep single ends if partner is discarded
	public static final String CFGKEY_THREADCOUNT = "threadcount"; //"auto" or number how many threads to use
	
	public static final String CFGKEY_IFILLUMINA = "ifillumina";
	public static final String CFGKEY_ISILLUMINA = "isillumina"; //--isIlluminaFormat=[<1.3/1.3/1.5/1.9/auto] ONLY IF ILLUMINA
	
	public static final String CFGKEY_CONVTOPHRED = "convtophred"; //--convertQualToPhred -> ONLY IF ILLUMINA
	public static final String CFGKEY_REMOVEN = "removen";
	
	public static final String CFGKEY_USEQUALTHRESHOLD = "usequalthreshold";//Average read-quality threshold
	public static final String CFGKEY_QUALTHRESHOLD = "qualthreshold";
	
	public static final String CFGKEY_USETRIMBYQUAL = "usetrimbyqual"; //Trim reads from left until threshold reached
	public static final String CFGKEY_TRIMBYQUAL = "trimbyqual";
	
	public static final String CFGKEY_TRIMBOTHENDS = "trimbothends";
	public static final String CFGKEY_OUTPUTFOLDER = "outputFolder";
	

	private final SettingsModelBoolean m_useotherfilterfile = new SettingsModelBoolean(CFGKEY_USEOTHERFILTERFILE, false);
	private final SettingsModelString m_otherfiltersettingsfile = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_OTHERFILTERSETTINGSFILE,"");

	private final SettingsModelBoolean m_removeadapters = new SettingsModelBoolean(CFGKEY_REMOVEADAPTERS, false);
	private final SettingsModelString m_adapters = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_ADAPTERS,"");
	
	private final SettingsModelBoolean m_dotrimpolyat = new SettingsModelBoolean(CFGKEY_DOTRIMPOLYAT, false);
	private final SettingsModelString m_trimpolyat = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_TRIMPOLYAT,"");
	
	private final SettingsModelBoolean m_filterfileexists = new SettingsModelBoolean(CFGKEY_FILTERFILEEXISTS, true);
	
	private final SettingsModelBoolean m_ifbarcodefile = new SettingsModelBoolean(CFGKEY_IFBARCODEFILE, false);
	private final SettingsModelString m_barcodefile = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_BARCODEFILE,"");
	
	private final SettingsModelBoolean m_lengthcutoff = new SettingsModelBoolean(CFGKEY_LENGTHCUTOFF, false);
	private final SettingsModelIntegerBounded m_minlength = new SettingsModelIntegerBounded(CFGKEY_MINLENGTH, 30, 1, Integer.MAX_VALUE);

	private final SettingsModelString m_preserve = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_PRESERVE,"");
	private final SettingsModelIntegerBounded m_threadcount = new SettingsModelIntegerBounded(CFGKEY_THREADCOUNT, 4, 1, Integer.MAX_VALUE);
	
	private final SettingsModelBoolean m_ifillumina = new SettingsModelBoolean(CFGKEY_IFILLUMINA, false);
	private final SettingsModelString m_isillumina = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_ISILLUMINA,"");
	
	private final SettingsModelString m_convtophred = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_CONVTOPHRED,"");
	private final SettingsModelString m_removen = new SettingsModelString(
			RawReadManipulatorNodeModel.CFGKEY_REMOVEN,"Remove reads containing Ns?");

	private final SettingsModelBoolean m_usequalthreshold = new SettingsModelBoolean(CFGKEY_USEQUALTHRESHOLD, false);
	private final SettingsModelIntegerBounded m_qualthreshold = new SettingsModelIntegerBounded(CFGKEY_QUALTHRESHOLD, 1, 1, Integer.MAX_VALUE);
	
	private final SettingsModelBoolean m_usetrimbyqual = new SettingsModelBoolean(CFGKEY_USETRIMBYQUAL, false);
	private final SettingsModelIntegerBounded m_trimbyqual = new SettingsModelIntegerBounded(CFGKEY_TRIMBYQUAL, 1, 1, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_trimbothends = new SettingsModelBoolean(CFGKEY_TRIMBOTHENDS, false);
	public final SettingsModelOptionalString m_outputfolder = new SettingsModelOptionalString(CFGKEY_OUTPUTFOLDER, "", false);
	
	//ReadType: paired-end or single-end
	private static String readType = "";
	
	
	/**
     * Constructor for the node model.
     */
    protected RawReadManipulatorNodeModel() {
    
        super(1, 1);

    	m_adapters.setEnabled(false);
		m_trimpolyat.setEnabled(false);
		m_minlength.setEnabled(false);
		m_preserve.setEnabled(true);
		m_isillumina.setEnabled(false);
		m_convtophred.setEnabled(false);
		m_barcodefile.setEnabled(false);
		m_qualthreshold.setEnabled(false);
		m_trimbyqual.setEnabled(false);
		m_otherfiltersettingsfile.setEnabled(false);
		m_useotherfilterfile.setEnabled(false);
		m_trimbothends.setEnabled(false);				
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	    	
    	String inFile1 = inData[0].iterator().next().getCell(0).toString();
    	String inFile2 = "";
    	
    	if(readType.equals("paired-end")){
    		inFile2 = inData[0].iterator().next().getCell(1).toString();
    	}
    	
    	
    	String outputFolder = m_outputfolder.getStringValue().isEmpty() ? "" : new File(m_outputfolder.getStringValue() + File.separator).getAbsolutePath();
  	
    	
    	/**Prepare Command**/
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("java");
    	command.add("-jar "+IO.getScriptPath()+"/libs/FastQReadFiltering.jar");

    	
    	String infileParameter = "--in="+inFile1;
    	if(readType.equals("paired-end") && !inFile2.equals("")) {
    		infileParameter+=","+inFile2; 
        	if(m_preserve.getStringValue().equals("Yes")) {
        		command.add("--preserve");
        	}
    	}
    	command.add(infileParameter);
    	
    	if(!outputFolder.isEmpty()) {
    		command.add("--outDir="+outputFolder);
    	}
    	
    	if(m_filterfileexists.getBooleanValue()){
    		command.add("--filtersettings="+inData[0].iterator().next().getCell(2).toString());
    	}
    	
		if(m_useotherfilterfile.getBooleanValue()) {
			command.add("--filtersettings="+m_otherfiltersettingsfile.getStringValue());
    	}
    	if(m_ifbarcodefile.getBooleanValue()){
    		command.add("--splitByAndTrimBarcodes="+m_barcodefile.getStringValue());
    	}
    	if(m_removeadapters.getBooleanValue()){
    		command.add("--removeAdapters="+m_adapters.getStringValue());
    	}
    	if(m_dotrimpolyat.getBooleanValue()){
    		String value = m_trimpolyat.getStringValue();
    		value = value.substring(0,1);
    		command.add("--trimPolyAT="+value);
    	}
    	if(m_lengthcutoff.getBooleanValue()){
    		command.add("--minLength="+m_minlength.getIntValue());
    	}
    	command.add("--threadCount="+m_threadcount.getIntValue());
    	if(m_ifillumina.getBooleanValue()){
    		command.add("--isIlluminaFormat="+m_isillumina.getStringValue());  
    		if(m_convtophred.getStringValue().equals("Yes")){
    			command.add("--convertQualToPhred");
    		}
    	}
    	if(m_removen.getStringValue().equals("Yes")){
    		command.add("--discardNReads");
    	}
    	if(m_usequalthreshold.getBooleanValue()){
    		command.add("--discardReadBelowAvgQuality="+m_qualthreshold.getIntValue());
    	}
    	if(m_usetrimbyqual.getBooleanValue()){
    		command.add("--trimByQual="+m_trimbyqual.getIntValue());
    	}
    	if(m_trimbothends.getBooleanValue()){
    		command.add("--trimByQualBothEnds");
    	}
    	
    	command.add("--noexit");
   	
    	String inFile1nozip = IO.removeZipExtension(inFile1);
    	String inFile2nozip = IO.removeZipExtension(inFile2);
    	
    	/** check if run was already sucessful **/
    	String[] com = command.toArray(new String[command.size()]);
    	File lockFile = new File(inFile1nozip.substring(0,inFile1nozip.lastIndexOf(".")) + ".RRM" +  SuccessfulRunChecker.LOCK_ENDING);

		String stdOutFile = inFile1nozip.substring(0,inFile1nozip.lastIndexOf(".")) + ".filtered.stdOut.log";
		String stdErrFile = inFile1nozip.substring(0,inFile1nozip.lastIndexOf(".")) + ".filtered.stdErr.log";
			
	    /**Execute**/
//	    Executor.executeCommand(new String[]{StringUtils.join(com, " ")}, exec, null, LOGGER, stdOutFile, stdErrFile, null);
	    super.executeCommand(new String[]{StringUtils.join(com, " ")}, exec, null, lockFile, stdOutFile, stdErrFile, null, null, null);

		
    	
        /**Create Output**/
	    
    	String outReadsFile1 = inFile1nozip.substring(0,inFile1nozip.lastIndexOf(".")) + ".filtered"+inFile1nozip.substring(inFile1nozip.lastIndexOf("."));
    	/**
         * fq format handling. fastq is already handled in Fastqc.jar
         */
        if (outReadsFile1.substring(outReadsFile1.length()-2, outReadsFile1.length()).equals("fq"))
    	{
        	outReadsFile1 = outReadsFile1.substring(0,inFile1nozip.lastIndexOf(".")) + ".fq.filtered.fastq";
    	}
        
    	if(!outputFolder.isEmpty()) outReadsFile1 = outputFolder + File.separator + new File(outReadsFile1).getName();
    	if(!new File(outReadsFile1).exists()) {	//If the file does not exist
			throw new Exception("The expected outfile "+outReadsFile1+" does not exist....something went wrong!");
    	}
    	
    	String outReadsFile2 = "";
    	if(readType.equals("paired-end")) {
    		outReadsFile2 = inFile2nozip.substring(0,inFile2nozip.lastIndexOf(".")) + ".filtered"+inFile2nozip.substring(inFile2nozip.lastIndexOf("."));
    		/**
             * fq format handling. fastq is already handled in Fastqc.jar
             */
            if (outReadsFile2.substring(outReadsFile2.length()-2, outReadsFile2.length()).equals("fq"))
        	{
            	outReadsFile2 = outReadsFile2.substring(0,inFile2nozip.lastIndexOf(".")) + ".fq.filtered.fastq";
        	}
            
    		if(!outputFolder.isEmpty()) outReadsFile1 = outputFolder + File.separator + new File(outReadsFile2).getName();
    		if(!new File(outReadsFile2).exists()) {
    			throw new IOException("The expected outfile "+outReadsFile2+" does not exist....something went wrong!");
        	}
    	}

    	BufferedDataContainer cont;
    	FileCell[] c;
    	
    	if(readType.equals("single-end")){
        	cont = exec.createDataContainer(createSpecs());
        	c = new FileCell[]{
        			(FileCell) FileCellFactory.create(outReadsFile1)};
    	}else{
        	cont = exec.createDataContainer(createSpecs());
        	c = new FileCell[]{
        			(FileCell) FileCellFactory.create(outReadsFile1),
        			(FileCell) FileCellFactory.create(outReadsFile2)};
    	}
    	
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
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	CompatibilityChecker CC = new CompatibilityChecker();
    	readType = CC.getReadType(inSpecs, 0);
    	if(CC.getWarningStatus()){
    		setWarningMessage(CC.getWarningMessages());
    	}
    	 	
    	// Check input ports
//    	String[] cn=inSpecs[0].getColumnNames();
//    	if(!cn[0].equals("") && !cn[0].equals("Path2ReadFile1")) {
//    		throw new InvalidSettingsException("This node is incompatible with the previous node. The outport of the previous node has to fit to the inport of this node.");
//    	}
    	    	
        return new DataTableSpec[]{createSpecs()};
    }

    
    /**
     * Create Tablespecs
     * @return
     */
    private DataTableSpec createSpecs(){
    	DataTableSpec out;
    	if(readType.equals("single-end")){ 	
    		out = new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL1, FastQCell.TYPE).createSpec()});
    	}else{
    		out = new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL1, FastQCell.TYPE).createSpec(),
        					new DataColumnSpecCreator(OUT_COL2, FastQCell.TYPE).createSpec()});
    	}
    	return out;
    }
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
    	super.saveSettingsTo(settings);
    	m_removeadapters.saveSettingsTo(settings);
    	m_adapters.saveSettingsTo(settings);
    	m_dotrimpolyat.saveSettingsTo(settings);
    	m_trimpolyat.saveSettingsTo(settings);
    	m_filterfileexists.saveSettingsTo(settings);
    	m_ifbarcodefile.saveSettingsTo(settings);
    	m_barcodefile.saveSettingsTo(settings);
    	m_lengthcutoff.saveSettingsTo(settings);
    	m_minlength.saveSettingsTo(settings);
    	m_preserve.saveSettingsTo(settings);
    	m_threadcount.saveSettingsTo(settings);
    	m_isillumina.saveSettingsTo(settings);
    	m_ifillumina.saveSettingsTo(settings);
    	m_convtophred.saveSettingsTo(settings);
    	m_removen.saveSettingsTo(settings);
    	m_usequalthreshold.saveSettingsTo(settings);
    	m_qualthreshold.saveSettingsTo(settings);
    	m_usetrimbyqual.saveSettingsTo(settings);
    	m_trimbyqual.saveSettingsTo(settings);
    	m_otherfiltersettingsfile.saveSettingsTo(settings);
    	m_useotherfilterfile.saveSettingsTo(settings);
    	m_trimbothends.saveSettingsTo(settings);
    	m_outputfolder.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    	super.loadValidatedSettingsFrom(settings);
    	m_removeadapters.loadSettingsFrom(settings);
    	m_adapters.loadSettingsFrom(settings);
    	m_dotrimpolyat.loadSettingsFrom(settings);
    	m_trimpolyat.loadSettingsFrom(settings);
    	m_filterfileexists.loadSettingsFrom(settings);
    	m_ifbarcodefile.loadSettingsFrom(settings);
    	m_barcodefile.loadSettingsFrom(settings);
    	m_lengthcutoff.loadSettingsFrom(settings);
    	m_minlength.loadSettingsFrom(settings);
    	m_preserve.loadSettingsFrom(settings);
    	m_threadcount.loadSettingsFrom(settings);
    	m_isillumina.loadSettingsFrom(settings);
    	m_ifillumina.loadSettingsFrom(settings);
    	m_convtophred.loadSettingsFrom(settings);
    	m_removen.loadSettingsFrom(settings);
    	m_usequalthreshold.loadSettingsFrom(settings);
    	m_qualthreshold.loadSettingsFrom(settings);
    	m_usetrimbyqual.loadSettingsFrom(settings);
    	m_trimbyqual.loadSettingsFrom(settings);
    	m_useotherfilterfile.loadSettingsFrom(settings);
    	m_otherfiltersettingsfile.loadSettingsFrom(settings);
    	m_trimbothends.loadSettingsFrom(settings);
    	m_outputfolder.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
   
    	super.validateSettings(settings);
    	m_removeadapters.validateSettings(settings);
    	m_adapters.validateSettings(settings);
    	m_dotrimpolyat.validateSettings(settings);
    	m_trimpolyat.validateSettings(settings);
    	m_filterfileexists.validateSettings(settings);
    	m_ifbarcodefile.validateSettings(settings);
    	m_barcodefile.validateSettings(settings);
    	m_lengthcutoff.validateSettings(settings);
    	m_minlength.validateSettings(settings);
    	m_preserve.validateSettings(settings);
    	m_threadcount.validateSettings(settings);
    	m_isillumina.validateSettings(settings);
    	m_ifillumina.validateSettings(settings);
    	m_convtophred.validateSettings(settings);
    	m_removen.validateSettings(settings);
    	m_usequalthreshold.validateSettings(settings);
    	m_qualthreshold.validateSettings(settings);
    	m_usetrimbyqual.validateSettings(settings);
    	m_trimbyqual.validateSettings(settings);
    	m_useotherfilterfile.validateSettings(settings);
    	m_otherfiltersettingsfile.validateSettings(settings);
    	m_trimbothends.validateSettings(settings);
    	m_outputfolder.validateSettings(settings);
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

}

