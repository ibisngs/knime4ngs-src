/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.vcfmerger;


import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
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
    		command.add("--variant "+IO.processFilePath(filename));
    	}
		
    	command.add("--genotypemergeoption "+m_GENOTYPEMERGEOPTION.getStringValue());
		
		return StringUtils.join(command, " ");
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	resetView();
    	this.FILES.clear();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	 super.saveSettingsTo(settings);
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
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
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
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
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
	protected String getOutfile() {
    	String Outfolder;
		try {
			Outfolder = IO.processFilePath(m_OUTFOLDER.getStringValue());
	    	String OUTFILETAG	= m_OUTFILETAG.getStringValue();
			String OUTFILE 		= Outfolder+"/"+OUTFILETAG+".vcf";
			return OUTFILE;
		} catch (InvalidSettingsException e) {
			e.printStackTrace();
		}

		return null;
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
    	for(String filename : FILES)
    	{
    		if (!filename.endsWith(".vcf")) {
    			throw new InvalidSettingsException("Inputfile is not a VCF file:"+filename);
    		}
    	}
	}

}

