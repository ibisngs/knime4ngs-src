package de.helmholtz_muenchen.ibis.ngs.mats;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Mats" Node.
 *
 * 
 * @author Michael Kluge
 */
public class MatsNodeDialog extends HTExecutorNodeDialog {

	private final static String BINARY_NAME = "MATS";
	
	
	
    /**
     * New pane for configuring the Mats node.
     */
    protected MatsNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString SET_OUTPUT_FOLDER			= new SettingsModelString(MatsNodeModel.CFGKEY_OUTPUT_FOLDER, MatsNodeModel.DEFAULT_OUTPUT_FOLDER);
    	final SettingsModelString SET_ANNOTATION_FILE		= new SettingsModelString(MatsNodeModel.CFGKEY_ANNOTATION_FILE, MatsNodeModel.DEFAULT_ANNOTATION_FILE);
    	final SettingsModelInteger SET_READ_LENGTH			= new SettingsModelInteger(MatsNodeModel.CFGKEY_READ_LENGTH, MatsNodeModel.DEFAULT_READ_LENGTH);	
    	final SettingsModelDouble SET_CUTOFF_DIFFERENCE		= new SettingsModelDouble(MatsNodeModel.CFGKEY_CUTOFF_DIFFERENCE, MatsNodeModel.DEFAULT_CUTOFF_DIFFERENCE);
    	final SettingsModelBoolean SET_ANALYSIS_TYPE		= new SettingsModelBoolean(MatsNodeModel.CFGKEY_ANALYSIS_TYPE, MatsNodeModel.DEFAULT_PAIRED_ANALYSIS);
    	final SettingsModelDouble SET_EXPRESSION_CHANGE		= new SettingsModelDouble(MatsNodeModel.CFGKEY_EXPRESSION_CHANGE, MatsNodeModel.DEFAULT_EXPRESSION_CHANGE);
        final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);

    	
    	DialogComponentFileChooser dcAnnotationFile = new DialogComponentFileChooser(SET_ANNOTATION_FILE, "his_id_INPUT_ANNOTATION_MATS", 0, ".gtf", ".GTF");
		DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER_Mats", 0, true);
		DialogComponentBoolean dcAnalysisType		= new DialogComponentBoolean(SET_ANALYSIS_TYPE, "paired analysis");
		DialogComponentNumber dcDiffCutoff			= new DialogComponentNumber(SET_CUTOFF_DIFFERENCE, "minimum difference splicing cutoff", 0.01);
		DialogComponentNumber dcReadLength			= new DialogComponentNumber(SET_READ_LENGTH, "read length", 1);
		DialogComponentNumber dcExpressionChange	= new DialogComponentNumber(SET_EXPRESSION_CHANGE, "gene expression foldchange filter", 0.5);
		
		// set a new title to them
		dcOutputFolder.setBorderTitle("Path to output folder");
		dcAnnotationFile.setBorderTitle("Path to gtf annotation file");
		
        addPrefPageSetting(SET_BINARY_PATH,getNameOfBinary());
		
		createNewTab("Mats Options");
		// add groups and components
		createNewGroup("Input");
		addDialogComponent(dcAnnotationFile);
		
		createNewGroup("Output");
		addDialogComponent(dcOutputFolder); 
		  
		createNewGroup("Further options");
		addDialogComponent(dcReadLength);
		addDialogComponent(dcAnalysisType);
		addDialogComponent(dcDiffCutoff);
		addDialogComponent(dcExpressionChange);
		
		// check range of diff
		SET_CUTOFF_DIFFERENCE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				// check range
				double value = SET_CUTOFF_DIFFERENCE.getDoubleValue();
				if(MatsNodeModel.MIN_CHANGE > value)
					SET_CUTOFF_DIFFERENCE.setDoubleValue(MatsNodeModel.MIN_CHANGE);
				if(MatsNodeModel.MAX_CHANGE <= value)
					SET_CUTOFF_DIFFERENCE.setDoubleValue(MatsNodeModel.MAX_CHANGE-0.001);
			}	
        });
		
		// check read length
		SET_READ_LENGTH.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				// check range
				if(SET_READ_LENGTH.getIntValue() <= 0)
					SET_READ_LENGTH.setIntValue(1);
			}	
        });
    }

	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}