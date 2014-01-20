package de.helmholtz_muenchen.ibis.ngs.samloader;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "SAMLoader" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class SAMLoaderNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the SAMLoader node.
     */
    protected SAMLoaderNodeDialog() {
    	
    	createNewGroup("Reference (e.g. genome) sequence: FastA file.");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(SAMLoaderNodeModel.CFGKEY_SEQFILE,null), "his0_id_samloader", 0, ""));
    	createNewGroup("Alignment/ SAM file.");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(SAMLoaderNodeModel.CFGKEY_SAMFILE,null), "his1_id_samloader", 0, "sam"));
    	
    }
}

