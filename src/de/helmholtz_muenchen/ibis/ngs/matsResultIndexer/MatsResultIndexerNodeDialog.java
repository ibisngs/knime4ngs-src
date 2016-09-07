package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "MatsResultIndexer" Node.
 * 
 * @author Michael Kluge
 */
public class MatsResultIndexerNodeDialog extends HTExecutorNodeDialog {

	private final static String BINARY_NAME = "index_gff";
    
    /**
     * New pane for configuring the MatsResultIndexer node.
     */
    protected MatsResultIndexerNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(MatsResultIndexerNodeModel.CFGKEY_OUTPUT_FILE, MatsResultIndexerNodeModel.DEFAULT_OUTPUT_FOLDER);
        final SettingsModelString SET_INPUT_FILE = new SettingsModelString(MatsResultIndexerNodeModel.CFGKEY_INPUT_FILE, MatsResultIndexerNodeModel.DEFAULT_INPUT_FOLDER);
        final SettingsModelBoolean SET_INCLUDE_NOVEL = new SettingsModelBoolean(MatsResultIndexerNodeModel.CFGKEY_INCLUDE_NOVEL, MatsResultIndexerNodeModel.DEFAULT_INCLUDE_NOVEL);
        final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);

		DialogComponentFileChooser dcInputFile 	= new DialogComponentFileChooser(SET_INPUT_FILE, "his_id_INPUT_FILE_MatsResultIndexer", 0, true);
		DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE_MatsResultIndexer", 0, true);
		DialogComponentBoolean dcIncludeNovel	 	= new DialogComponentBoolean(SET_INCLUDE_NOVEL, "include novel events found by MATS");
		  	
		// set a new title to them
		dcOutputFile.setBorderTitle("Path to output folder");
		dcInputFile.setBorderTitle("Path to input folder");
		
        addPrefPageSetting(SET_BINARY_PATH,getNameOfBinary());
		
		createNewTab("MatsResultIndexer");
		    // add groups and components
		createNewGroup("Input");
		addDialogComponent(dcInputFile);
		 
		createNewGroup("Output");
		addDialogComponent(dcOutputFile);
		 
		createNewGroup("Further options");
		addDialogComponent(dcIncludeNovel);
    }
    
	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}

