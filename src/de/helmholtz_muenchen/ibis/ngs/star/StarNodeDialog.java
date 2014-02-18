package de.helmholtz_muenchen.ibis.ngs.star;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeDialog extends DefaultNodeSettingsPane {

	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH 		= StarNodeModel.getSettingsModelString(StarNodeModel.CFGKEY_BINARY_PATH, this);
    private final SettingsModelString SET_PARAMETER_FILE 	= StarNodeModel.getSettingsModelString(StarNodeModel.CFGKEY_PARAMETER_FILE, this);
    private final SettingsModelString SET_RUN_MODE 			= StarNodeModel.getSettingsModelString(StarNodeModel.CFGKEY_RUN_MODE, this);
    private final SettingsModelString SET_OUTPUT_FOLDER 	= StarNodeModel.getSettingsModelString(StarNodeModel.CFGKEY_OUTPUT_FOLDER, this);
    private final SettingsModelString SET_GENOME_FOLDER 	= StarNodeModel.getSettingsModelString(StarNodeModel.CFGKEY_GENOME_FOLDER, this);
    
    protected StarNodeDialog() {
        super();
       
        // create open file/folder components
        DialogComponentFileChooser dcBinaryPath 	= new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BINARY_PATH", 0);
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER", 0, true);
        DialogComponentFileChooser dcParameterFile 	= new DialogComponentFileChooser(SET_PARAMETER_FILE, "his_id_PARAMETER_FILE", 0);
       	DialogComponentFileChooser dcGenomeFolder 	= new DialogComponentFileChooser(SET_GENOME_FOLDER, "his_id_GENOME_FOLDER", 0, true);

       	// create string selection component
       	DialogComponentStringSelection dcRunMode 	= new DialogComponentStringSelection(SET_RUN_MODE, "runMode:", StarNodeModel.DEFAULT_RUN_MODE, StarNodeModel.ALTERNATIVE_RUN_MODE);
       	
       	// set a new title to them
       	dcBinaryPath.setBorderTitle("path to STAR binary");
       	dcOutputFolder.setBorderTitle("path to output folder");
       	dcParameterFile.setBorderTitle("path to STAR parameter file");
       	dcGenomeFolder.setBorderTitle("path to indexed genome");
     
       	// add groups and components
        createNewGroup("generell options");
        addDialogComponent(dcRunMode);
        addDialogComponent(dcBinaryPath);
        
        createNewGroup("input");
        addDialogComponent(dcParameterFile);
        addDialogComponent(dcGenomeFolder);
        
        createNewGroup("output");
        addDialogComponent(dcOutputFolder);
 
        // add change listener to runMode
        SET_RUN_MODE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(StarNodeModel.DEFAULT_RUN_MODE.equals(SET_RUN_MODE.getStringValue()))
					SET_GENOME_FOLDER.setEnabled(true);
				else
					SET_GENOME_FOLDER.setEnabled(false);
			}
        });
    }
}

