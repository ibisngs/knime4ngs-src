package de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants;


import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

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
public abstract class SelectVariantsNodeDialog extends GATKNodeDialog {

	/**
	 * The SettingsModels
	 */
	
//    private final SettingsModelString m_FILTERSTRING; 
    
    /**
     * New pane for configuring the GATKSelectVariants node.
     */
//    protected SelectVariantsNodeDialog() {
//
//    	createNewGroup("Path to GATK jar file");
//    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
//    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
//    	addDialogComponent(gatkf);
//    	
//    	createNewGroup("Reference Genome");
//    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
//    	gatkf.setBorderTitle("Choose the reference genome");
//    	addDialogComponent(ref_genome);
//    	
//        createNewGroup("");
//        addDialogComponent(new DialogComponentNumber(VCFCOLUMN, "Select Column of inData in which Path to Input VCF is specified", 1));
//    	
//        addFilterDialogComponent();
//        
//    }
    
    protected void addDialogComponent() {
    	
        SettingsModelIntegerBounded VCFCOLUMN = new SettingsModelIntegerBounded(SelectVariantsNodeModel.CFGKEY_VCFCOLUMN, 1, 1, Integer.MAX_VALUE);
    	createNewTab("SelectVariants");
        createNewGroup("");
        addDialogComponent(new DialogComponentNumber(VCFCOLUMN, "Column of input VCF file", 1));
        addFilterDialogComponent();
    }
    
    protected SettingsModelString getFILTERSTRINGModel(){
    	return new SettingsModelString(SelectVariantsNodeModel.CFGKEY_FILTERSTRING, "");
    }
    
    /**
     * Implement specific filter gui
     */
    protected abstract void addFilterDialogComponent();
    
}

