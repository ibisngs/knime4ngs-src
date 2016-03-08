package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import java.io.File;
import java.util.LinkedHashMap;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;

/**
 * This is the model implementation of MatsResultIndexer.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultIndexerNodeModel extends BinaryWrapperNodeModel {
	
	public static final String OUT_COL1 = "OutputFolder";
	public static final String OUT_COL2 = "OutputGFF3File";
	public static final String GFF_FILE_NAME = "events4miso.gff3";
	
	 // keys for SettingsModels
    protected static final String CFGKEY_OUTPUT_FILE 			= "OutputFolder";
    protected static final String CFGKEY_INPUT_FILE				= "InputFolder";
    protected static final String CFGKEY_INCLUDE_NOVEL			= "IncludeNovel";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_OUTPUT_FOLDER 		= "./output/";	
    protected static final String DEFAULT_INPUT_FOLDER 			= ".";
    protected static final boolean DEFAULT_INCLUDE_NOVEL		= true;
    
    private final static String NAME_OF_INPUT_FILE 				= "--index";
    private final static String NAME_OF_OUTPUT_FILE 			= " ";	
    
    private final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(CFGKEY_OUTPUT_FILE, DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelString SET_INPUT_FILE = new SettingsModelString(CFGKEY_INPUT_FILE, DEFAULT_INPUT_FOLDER);
    private final SettingsModelBoolean SET_INCLUDE_NOVEL = new SettingsModelBoolean(CFGKEY_INCLUDE_NOVEL, DEFAULT_INCLUDE_NOVEL);
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(MatsResultIndexerNodeModel.class);
	
	 protected MatsResultIndexerNodeModel() {
        super(0, 1, false, false);
        this.addSetting(SET_INPUT_FILE);
    	this.addSetting(SET_OUTPUT_FILE);
    	this.addSetting(SET_INCLUDE_NOVEL);
    }

	@Override
	protected LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData) {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();
		
		/********************* SIMPLE PARAMETER ***************************/		
		pars.put(NAME_OF_INPUT_FILE, SET_OUTPUT_FILE.getStringValue() + File.separator + GFF_FILE_NAME);
		pars.put(NAME_OF_OUTPUT_FILE, SET_OUTPUT_FILE.getStringValue());
		
		return pars;
	}
	
	 protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		 // create index GFF file first
		 String inputFolder = SET_INPUT_FILE.getStringValue();
		 String outputFile = SET_OUTPUT_FILE.getStringValue();
		 
		 // create folder, if not existent
		 File out = new File(outputFile);
		 outputFile = outputFile + File.separator + GFF_FILE_NAME;
		 if(!out.exists())
			 out.mkdirs();
		 
		 // parse the annotation
		 if(new File(inputFolder).exists() & new File(inputFolder).isDirectory()) {
			 @SuppressWarnings("unused")
			 IndexCreator ic = new IndexCreator(inputFolder, SET_INCLUDE_NOVEL.getBooleanValue(), outputFile);
			 
			 // check, if output file is there
			 if(!new File(outputFile).exists()) {
				 throw new IllegalArgumentException("Parsing of '"+ inputFolder +"' failed. Make sure that it is a valid ASEvents folder of MATS.");
			 }
			 
			 return super.execute(inData, exec);
		 }
		 else {
			 throw new IllegalArgumentException("Folder '"+ inputFolder +"' was not found.");
		 }
	 }
	
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, StringCell.TYPE).createSpec()});
    }

	@Override
	protected BufferedDataTable[] getOutputData(ExecutionContext exec, String command, BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(getDataOutSpec1());
		
    	DataCell[] c = new DataCell[]{
    			new StringCell(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true)),
    			new StringCell(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true) + File.separator + GFF_FILE_NAME)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
		
        return new BufferedDataTable[]{cont.getTable()};
	}

	@Override
	protected File getPathToStderrFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true) + "MatsResultIndexer.out");
	}

	@Override
	protected File getPathToStdoutFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true) + "MatsResultIndexer.err");
	}

	@Override
	protected File getPathToLockFile() {
		return null;
	}
}

