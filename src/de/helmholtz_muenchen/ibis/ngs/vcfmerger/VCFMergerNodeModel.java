package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;

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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;


/**
 * This is the model implementation of VCFMerger.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeModel extends NodeModel {
    
    public static final String CFGKEY_INFOLDER					= "infolder";
    public static final String CFGKEY_REGEX						= "regex";
    public static final String CFGKEY_OUTFOLDER					= "outfolder";
    public static final String CFGKEY_GATK						= "gatk";
	public static final String CFGKEY_REF_GENOME 				= "REFGENOME";
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
	public static final String CFGKEY_OUTFILETAG		 		= "outfiletag";
    
    private final SettingsModelString m_GATK 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REF_GENOME, "");
//    private final SettingsModelString m_INFOLDER 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_INFOLDER, "");
//    private final SettingsModelString m_REGEX 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REGEX, "");
    private final SettingsModelString m_OUTFOLDER 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");
    private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
    private final SettingsModelString m_OUTFILETAG			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFILETAG, "");

    
    // keys for SettingsModels
    protected static final String CFGKEY_FILE_LIST 			= "FileList";
    protected static final String CFGKEY_FILE_DIR 			= "FileDir";
    protected static final String CFGKEY_FILE_FILE 			= "FileFile";
    protected static final String CFGKEY_FILE_LIST_DISPLAY 	= "FileListDisplay";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_FILE_DIR = "-";
    protected static final String DEFAULT_FILE_FILE = "-";
    protected static final String DEFAULT_REGEX = ".*";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString m_SET_FILE_DIR 				= new SettingsModelString(CFGKEY_FILE_DIR, DEFAULT_FILE_DIR);
    private final SettingsModelString m_SET_FILE_FILE 				= new SettingsModelString(CFGKEY_FILE_FILE, DEFAULT_FILE_FILE);
    private final SettingsModelString m_SET_NAME_REGEX				= new SettingsModelString(CFGKEY_REGEX, DEFAULT_REGEX);
    
    // storage 
    private final HashSet<String> FILES		= new HashSet<String>();
    private boolean hasConfigureOpendOnce 	= false; // true, if configure was opend once
    
	public static final String FILENAME_END_REGEX 		= "\\.(vcf)$";
    
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFMergerNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1 = "MERGED VARIANTS";
	
    /**
     * Constructor for the node model.
     */
    protected VCFMergerNodeModel() {
    
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

    	String Outfolder 	= m_OUTFOLDER.getStringValue();
    	String OUTFILETAG	= m_OUTFILETAG.getStringValue();
    	
		String OUTFILE = Outfolder+"/AllSamples_vcfmerger"+OUTFILETAG+".vcf";
		String ERRFILE = Outfolder+"/AllSamples_vcfmerger"+OUTFILETAG+".vcf.err";
		
		command.add("java");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T CombineVariants");
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	for(String filename : FILES)
    	{
    		command.add("--variant "+filename);
    	}
		command.add("--genotypemergeoption "+m_GENOTYPEMERGEOPTION.getStringValue());
		command.add("-o "+OUTFILE);
		
		try {
			Executor.executeCommand(new String[]{StringUtils.join(command, " ")}, exec,new String[]{}, LOGGER,OUTFILE,ERRFILE);
			
		} catch (CanceledExecutionException | InterruptedException
				| ExecutionException | UnsuccessfulExecutionException e) {
			e.printStackTrace();
		}
    	    	
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
    	this.FILES.clear();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	if(hasConfigureOpendOnce && FILES.size() == 0)
    		throw new InvalidSettingsException("Select at least one VCF file.");
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
//         m_INFOLDER.saveSettingsTo(settings);
         m_OUTFOLDER.saveSettingsTo(settings);
         m_GATK.saveSettingsTo(settings);
         m_REF_GENOME.saveSettingsTo(settings);
         m_GENOTYPEMERGEOPTION.saveSettingsTo(settings);
         m_OUTFILETAG.saveSettingsTo(settings);
         
         m_SET_FILE_DIR.saveSettingsTo(settings);
         m_SET_FILE_FILE.saveSettingsTo(settings);
         m_SET_NAME_REGEX.saveSettingsTo(settings);
        
     	settings.addStringArray(VCFMergerNodeModel.CFGKEY_FILE_LIST, FILES.toArray(new String[FILES.size()]));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
//        m_INFOLDER.loadSettingsFrom(settings);
        m_OUTFOLDER.loadSettingsFrom(settings);
//        m_REGEX.loadSettingsFrom(settings);
        m_GATK.loadSettingsFrom(settings);
        m_REF_GENOME.loadSettingsFrom(settings);
        m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
        m_OUTFILETAG.loadSettingsFrom(settings);
        
        m_SET_FILE_DIR.loadSettingsFrom(settings);
        m_SET_FILE_FILE.loadSettingsFrom(settings);
        m_SET_NAME_REGEX.loadSettingsFrom(settings);
        
    	// clean the old data
    	FILES.clear();
    	// check, if data is set
        if (settings.containsKey(VCFMergerNodeModel.CFGKEY_FILE_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(VCFMergerNodeModel.CFGKEY_FILE_LIST))
					this.FILES.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
//        m_INFOLDER.validateSettings(settings);
        m_OUTFOLDER.validateSettings(settings);
//        m_REGEX.validateSettings(settings);
        m_GATK.validateSettings(settings);
        m_REF_GENOME.validateSettings(settings);
        m_GENOTYPEMERGEOPTION.validateSettings(settings);
        m_OUTFILETAG.validateSettings(settings);
        
        m_SET_FILE_DIR.validateSettings(settings);
        m_SET_FILE_FILE.validateSettings(settings);
        m_SET_NAME_REGEX.validateSettings(settings);
        
        // configure must have been opened or we won't be here
        hasConfigureOpendOnce = true;
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

