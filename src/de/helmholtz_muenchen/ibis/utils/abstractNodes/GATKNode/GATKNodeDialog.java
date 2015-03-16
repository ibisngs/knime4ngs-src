package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import javax.swing.JFileChooser;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponent;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObjectSpec;


public abstract class GATKNodeDialog extends DefaultNodeSettingsPane{

	
	
    private final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);

	
    protected GATKNodeDialog() {
    	
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	addDialogComponent(new DialogComponentNumber(m_GATK_MEM, "GATK Memory", 1));
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	ref_genome.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
    	addDialogComponent();

    }
	
	
    protected abstract void addDialogComponent();
    
}
