package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulator;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.fastqc.FastQCNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;


/**
 * This is the model implementation of RawReadManipulator.
 * 
 *
 * @author Sebastian
 */
public class RawReadManipulatorNodeModel extends NodeModel {
    
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
    	String inFile2 = inData[0].iterator().next().getCell(1).toString();
    	//String filterFile = inData[0].iterator().next().getCell(2).toString();
    	String readType = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String path = FastQCNodeModel.class.getProtectionDomain().getCodeSource().getLocation().getPath();
    	
    	/**Initialize logfile**/
    	String logfile = inFile1.substring(0,inFile1.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(ShowOutput.getNodeStartTime("RawReadManipulator"));
    	ShowOutput.writeLogFile(logBuffer);
    	/**end initializing logfile**/
    	
       	/**Construct string with all parameters**/
    	StringBuffer call = new StringBuffer(32);
    	call.append("--in="+inFile1);
    	if(readType.equals("paired-end") && !inFile2.equals("")) {
    		call.append(","+inFile2);   
        	if(m_preserve.getStringValue().equals("Yes")) {
        		call.append(" --preserve");
        	}
    	}
    	if(m_filterfileexists.getBooleanValue()){
    		call.append(" --filtersettings="+inData[0].iterator().next().getCell(2).toString());
    	}
		if(m_useotherfilterfile.getBooleanValue()) {
			call.append(" --filtersettings="+m_otherfiltersettingsfile.getStringValue());
    	}
    	if(m_ifbarcodefile.getBooleanValue()){
    		call.append(" --splitByAndTrimBarcodes="+m_barcodefile.getStringValue());
    	}
    	if(m_removeadapters.getBooleanValue()){
    		call.append(" --removeAdapters="+m_adapters.getStringValue());
    	}
    	if(m_dotrimpolyat.getBooleanValue()){
    		String value = m_trimpolyat.getStringValue();
    		value = value.substring(0,1);
    		call.append(" --trimPolyAT="+value);
    	}
    	if(m_lengthcutoff.getBooleanValue()){
    		call.append(" --minLength="+m_minlength.getIntValue());
    	}
    	call.append(" --threadCount="+m_threadcount.getIntValue());
    	if(m_ifillumina.getBooleanValue()){
    		call.append(" --isIlluminaFormat="+m_isillumina.getStringValue());  
    		if(m_convtophred.getStringValue().equals("Yes")){
    			call.append(" --convertQualToPhred");
    		}
    	}
    	if(m_removen.getStringValue().equals("Yes")){
    		call.append(" --discardNReads");
    	}
    	if(m_usequalthreshold.getBooleanValue()){
    		call.append(" --discardReadBelowAvgQuality="+m_qualthreshold.getIntValue());
    	}
    	if(m_usetrimbyqual.getBooleanValue()){
    		call.append(" --trimByQual="+m_trimbyqual.getIntValue());
    	}
    	if(m_trimbothends.getBooleanValue()){
    		call.append(" --trimByQualBothEnds");
    	}
    	
    	call.append(" --noexit");
    	//String[] callReady = call.toString().split(" ");
   	
    	FileOutputStream errFile = new FileOutputStream(new File(ShowOutput.getLogFile()), true);
    	FileOutputStream outFile = new FileOutputStream(new File(ShowOutput.getLogFile()), true);
    	
    	PrintStream stdErr = new PrintStream(errFile); //append
    	PrintStream stdOut = new PrintStream(outFile);
    	System.setOut(stdOut);
    	System.setErr(stdErr);
    	
    	String com = "java -jar "+path+"/libs/FastQReadFiltering.jar "+call;
    	System.out.println(com);
    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
    	Process p1 = b.start();
    	p1.waitFor();
        logBuffer.append(ShowOutput.getLogEntry(p1, com));
        
    	//RawReadManipulator.main(callReady);

    	/*if(readType.equals("paired-end") && !inFile2.equals("")) {
    		callReady[0] = "--in="+inFile2;
   		if(callReady[1].lastIndexOf("filtersettings") != -1) {
    			callReady[1] = "--filtersettings="+inData[0].iterator().next().getCell(3).toString();
    		}
    		RawReadManipulator.main(callReady);
    	}*/
    	
        //Create Output
    	String outReadsFile1 = inFile1 + ".filtered.fastq";
    	if(!new File(outReadsFile1).exists()) {
    		outReadsFile1 = inFile1.substring(0,inFile1.lastIndexOf(".")) + ".filtered.fastq";
    	}
    	String outReadsFile2 = "";
    	if(!inFile2.equals("")) {
    		outReadsFile2 = inFile2 + ".filtered.fastq";
    		if(!new File(outReadsFile2).exists()) {
        		outReadsFile2 = inFile2.substring(0,inFile2.lastIndexOf(".")) + ".filtered.fastq";
        	}
    	}
        DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2ReadFile1", StringCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2ReadFile2", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(outReadsFile1);
    	StringCell cl2 = new StringCell(outReadsFile2);
    	DataCell[] c = new DataCell[]{cl1,cl2};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();

    	outFile.flush();
    	errFile.flush();
    	stdOut.close();
    	stdErr.close();
    	ShowOutput.writeLogFile(new StringBuffer(ShowOutput.getNodeEndTime()));
    	
        return new BufferedDataTable[]{out};
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
    	
    	
    	// Check input ports
    	String[] cn=inSpecs[0].getColumnNames();
    	if(!cn[0].equals("") && !cn[0].equals("Path2ReadFile1")) {
    		throw new InvalidSettingsException("This node is incompatible with the previous node. The outport of the previous node has to fit to the inport of this node.");
    	}
    	
    	/*//Check if output file exists
    	String outfile = "";
    	File outpath = new File(outfile);
        if(outpath.exists()){
        	setWarningMessage(outfile+" already exists! Please rename or move to other directory.");
        }*/
    	
        return new DataTableSpec[]{null};
    }

    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
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

