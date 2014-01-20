package de.helmholtz_muenchen.ibis.ngs.bamloader;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "BAMLoader" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class BAMLoaderNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the BAMLoader node.
     */
    protected BAMLoaderNodeDialog() {
    	
    	createNewGroup("SamTools.");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(BAMLoaderNodeModel.CFGKEY_SAMTOOLS,null), "his_baml_id", 0, ""));
    	createNewGroup("Reference (e.g. genome) sequence: FastA file.");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(BAMLoaderNodeModel.CFGKEY_SEQFILE,null), "his0_baml_id", 0, ""));
    	createNewGroup("Alignment/ BAM file (sorted, indexed, binary).");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(BAMLoaderNodeModel.CFGKEY_BAMFILE,null), "his1_baml_id", 0, "bam"));
    	
    }
}

