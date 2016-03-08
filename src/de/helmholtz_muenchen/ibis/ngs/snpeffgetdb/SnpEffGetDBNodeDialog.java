package de.helmholtz_muenchen.ibis.ngs.snpeffgetdb;

import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SnpEffGetDB" Node.
 * 
 * @author Maximilian Hastreiter 
 */
public class SnpEffGetDBNodeDialog extends HTExecutorNodeDialog {

	
	final SettingsModelString snpeff_folder = new SettingsModelString(
			SnpEffGetDBNodeModel.CFGKEY_SNPEFF_FOLDER, null);
	final SettingsModelString database = new SettingsModelString(
			SnpEffGetDBNodeModel.CFGKEY_DATABASE, null);
	//final SettingsModelString database_folder = new SettingsModelString(
	//		SnpEffGetDBNodeModel.CFGKEY_DATABASE_FOLDER, null);
    /**
     * New pane for configuring the SnpEffGetDB node.
     */
    protected SnpEffGetDBNodeDialog() {
    	
    	super();
    	
    	createNewGroup("snpEff directory");
    	addDialogComponent(new DialogComponentFileChooser(snpeff_folder, "par_1", 0, true));
    	
    	createNewGroup("Database name");
    	addDialogComponent(new DialogComponentString(database, ""));
    	
    	//createNewGroup("Database folder (where to extract the database)");
    	//addDialogComponent(new DialogComponentFileChooser(database_folder, "par_2", 0, true));
    	
    }
}

