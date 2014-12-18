package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


public abstract class GATKNodeDialog extends DefaultNodeSettingsPane{

	
	
    private final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
	
	
    protected GATKNodeDialog() {
    	
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	gatkf.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
    	addDialogComponent();
    	
    }
	
	
    protected abstract void addDialogComponent();
    
}
