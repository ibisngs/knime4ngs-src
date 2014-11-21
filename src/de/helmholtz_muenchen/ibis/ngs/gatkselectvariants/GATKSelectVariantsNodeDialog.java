package de.helmholtz_muenchen.ibis.ngs.gatkselectvariants;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKSelectVariants" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class GATKSelectVariantsNodeDialog extends DefaultNodeSettingsPane {

	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString GATK = new SettingsModelString(GATKSelectVariantsNodeModel.CFGKEY_GATK_PATH, "");
    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKSelectVariantsNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded VCFCOLUMN = new SettingsModelIntegerBounded(GATKSelectVariantsNodeModel.CFGKEY_VCFCOLUMN, 1, 1, Integer.MAX_VALUE);
	
	
    /**
     * New pane for configuring the GATKSelectVariants node.
     */
    protected GATKSelectVariantsNodeDialog() {

    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	gatkf.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
        createNewGroup("");
        addDialogComponent(new DialogComponentNumber(VCFCOLUMN, "Select Column of inData in which Path to Input VCF is specified", 1));
    	
    }
}

