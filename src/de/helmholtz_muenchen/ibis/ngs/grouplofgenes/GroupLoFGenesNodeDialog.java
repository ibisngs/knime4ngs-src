package de.helmholtz_muenchen.ibis.ngs.grouplofgenes;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GroupLoFGenes" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class GroupLoFGenesNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString gene_summary = new SettingsModelString(GroupLoFGenesNodeModel.CFGKEY_GENESUM,"-");
	private final SettingsModelString cds_file = new SettingsModelString(GroupLoFGenesNodeModel.CFGKEY_CDSFILE, "-");
	private final SettingsModelInteger samples = new SettingsModelInteger(GroupLoFGenesNodeModel.CFGKEY_SAMPLES,0);
	
    /**
     * New pane for configuring the GroupLoFGenes node.
     */
    protected GroupLoFGenesNodeDialog() {
    	
    	createNewGroup("Select gene summary");
    	addDialogComponent(new DialogComponentFileChooser(gene_summary, "his_id_gene_summary", 0, ".gene_summary.tsv"));

    	createNewGroup("Select CDS file");
    	addDialogComponent(new DialogComponentFileChooser(cds_file, "his_id_cds_file", 0, ".fa|.fasta"));
    	
    	createNewGroup("Study details");
    	addDialogComponent(new DialogComponentNumber(samples, "Number of samples",0));
    }
}

