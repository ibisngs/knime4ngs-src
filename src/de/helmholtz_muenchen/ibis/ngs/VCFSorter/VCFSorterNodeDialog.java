package de.helmholtz_muenchen.ibis.ngs.VCFSorter;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "VCFSorter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Kaarin Ahomaa
 */
public class VCFSorterNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VCFSorter node.
     */
	private final SettingsModelString refseq = new SettingsModelString(VCFSorterNodeModel.CFGKEY_REFSEQFILE,"");
	
	
    protected VCFSorterNodeDialog() {
    	
    	createNewGroup("Reference sequence: FastA file (e.g. genome)");
    	addDialogComponent(new DialogComponentFileChooser(refseq, "his1_id_BWA", 0, ".fa|.fasta"));
    	
    	
    }
}

