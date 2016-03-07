package de.helmholtz_muenchen.ibis.ngs.star;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeDialog;

/**
 * <code>NodeDialog</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeDialog extends BinaryWrapperNodeDialog {

	private final static String BINARY_NAME = "STAR";
	
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_RUN_MODE			= new SettingsModelString(StarNodeModel.CFGKEY_RUN_MODE, StarNodeModel.DEFAULT_RUN_MODE);
    private final SettingsModelString SET_OUTPUT_FOLDER		= new SettingsModelString(StarNodeModel.CFGKEY_OUTPUT_FOLDER, StarNodeModel.DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelString SET_GENOME_FOLDER		= new SettingsModelString(StarNodeModel.CFGKEY_GENOME_FOLDER, StarNodeModel.DEFAULT_GENOME_FOLDER);
    
    protected StarNodeDialog() {
        super();
       
        // create open file/folder components
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER", 0, true);
       	DialogComponentFileChooser dcGenomeFolder 	= new DialogComponentFileChooser(SET_GENOME_FOLDER, "his_id_GENOME_FOLDER", 0, true);

       	// create string selection component
       	DialogComponentStringSelection dcRunMode 	= new DialogComponentStringSelection(SET_RUN_MODE, "runMode:", StarNodeModel.DEFAULT_RUN_MODE, StarNodeModel.ALTERNATIVE_RUN_MODE);
       	
       	// set a new title to them
       	dcOutputFolder.setBorderTitle("Path to output folder");
       	dcGenomeFolder.setBorderTitle("Path to indexed genome");
     
       	// add groups and components
       	createNewGroup("STAR Options");
        addDialogComponent(dcRunMode);

        createNewGroup("Input");
        addDialogComponent(dcGenomeFolder);
        
        createNewGroup("Output");
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

	@Override
	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}

