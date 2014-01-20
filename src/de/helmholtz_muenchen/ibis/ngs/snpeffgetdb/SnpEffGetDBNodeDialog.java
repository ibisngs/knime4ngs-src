package de.helmholtz_muenchen.ibis.ngs.snpeffgetdb;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "SnpEffGetDB" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class SnpEffGetDBNodeDialog extends DefaultNodeSettingsPane {

	
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

    	createNewGroup("snpEff directory");
    	addDialogComponent(new DialogComponentFileChooser(snpeff_folder, "par_1", 0, true));
    	
    	createNewGroup("Database name");
    	addDialogComponent(new DialogComponentString(database, ""));
    	
    	//createNewGroup("Database folder (where to extract the database)");
    	//addDialogComponent(new DialogComponentFileChooser(database_folder, "par_2", 0, true));
    	
    }
}

