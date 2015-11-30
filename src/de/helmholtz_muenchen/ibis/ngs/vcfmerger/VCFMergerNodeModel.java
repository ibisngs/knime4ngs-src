package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;


/**
 * This is the model implementation of VCFMerger.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeModel extends GATKNodeModel {
    
    public static final String CFGKEY_INFOLDER					= "infolder";
    public static final String CFGKEY_REGEX						= "regex";
    public static final String CFGKEY_OUTFOLDER					= "outfolder";
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
	public static final String CFGKEY_OUTFILETAG		 		= "outfiletag";
    
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
    protected VCFMergerNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
   	}

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		ArrayList<String> command 	= new ArrayList<String>();
		
    	for(String filename : FILES)
    	{
    		command.add("--variant "+filename);
    	}
		
    	command.add("--genotypemergeoption "+m_GENOTYPEMERGEOPTION.getStringValue());
		
		return StringUtils.join(command, " ");
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
    protected void saveExtraSettingsTo(final NodeSettingsWO settings) {
         m_OUTFOLDER.saveSettingsTo(settings);
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
    protected void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_OUTFOLDER.loadSettingsFrom(settings);
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
    protected void validateExtraSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_OUTFOLDER.validateSettings(settings);
        m_GENOTYPEMERGEOPTION.validateSettings(settings);
        m_OUTFILETAG.validateSettings(settings);     
        m_SET_FILE_DIR.validateSettings(settings);
        m_SET_FILE_FILE.validateSettings(settings);
        m_SET_NAME_REGEX.validateSettings(settings);
        
        // configure must have been opened or we won't be here
        hasConfigureOpendOnce = true;
    }
    
	@Override
	protected String getCommandWalker() {
		return "CombineVariants";
	}

	@Override
	protected File getLockFile() {
		return new File(getOutfile()+SuccessfulRunChecker.LOCK_ENDING);
	}

	@Override
	protected String getOutfile() {
    	String Outfolder 	= m_OUTFOLDER.getStringValue();
    	String OUTFILETAG	= m_OUTFILETAG.getStringValue();
		String OUTFILE 		= Outfolder+"/"+OUTFILETAG+".vcf";
		return OUTFILE;
	}

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		return true;
	}

	@Override
	protected DataType getOutColType() {
		return VCFCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
    	if(hasConfigureOpendOnce && FILES.size() == 0)
    		throw new InvalidSettingsException("Select at least one VCF file.");
		
	}

}

