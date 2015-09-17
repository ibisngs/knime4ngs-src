package de.helmholtz_muenchen.ibis.ngs.gatkcombinegvcfs;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "CombineGVCFs" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class CombineGVCFsNodeDialog extends GATKNodeDialog {

	@Override
	protected void addDialogComponent() {
		
		final SettingsModelString BED_FILE 				= new SettingsModelString(CombineGVCFsNodeModel.CFGKEY_BED_FILE, "");
	    final SettingsModelBoolean BED_FILE_CHECKBOX 	= new SettingsModelBoolean(CombineGVCFsNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
			    
    	createNewGroup("BED File");
    	addDialogComponent(new DialogComponentBoolean(BED_FILE_CHECKBOX, "Use BED file?"));
    	addDialogComponent(new DialogComponentFileChooser(BED_FILE, "BED_FILE", JFileChooser.OPEN_DIALOG, false, ".bed"));
		
		BED_FILE_CHECKBOX.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					BED_FILE.setEnabled(BED_FILE_CHECKBOX.getBooleanValue());
				}
		});
	}
}

