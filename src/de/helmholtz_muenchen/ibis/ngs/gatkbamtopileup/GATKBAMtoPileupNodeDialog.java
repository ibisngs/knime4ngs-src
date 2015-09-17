package de.helmholtz_muenchen.ibis.ngs.gatkbamtopileup;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKBAMtoPileup" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKBAMtoPileupNodeDialog extends DefaultNodeSettingsPane {

	 	private final SettingsModelString GATK 				= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_GATK_PATH, "---");
//	    private final SettingsModelString INFILE 			= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_INFILE, "---");
	    private final SettingsModelString REF_GENOME 		= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_REF_GENOME, "---");
//	    private final SettingsModelOptionalString BED_FILE	= new SettingsModelOptionalString(GATKBAMtoPileupNodeModel.CFGKEY_BED_FILE, "---",false);
	    private final SettingsModelIntegerBounded GATK_NT 	= new SettingsModelIntegerBounded(GATKBAMtoPileupNodeModel.CFGKEY_GATK_NT, 1, 1, Integer.MAX_VALUE);
	    private final SettingsModelIntegerBounded GATK_NCT 	= new SettingsModelIntegerBounded(GATKBAMtoPileupNodeModel.CFGKEY_GATK_NCT, 1, 1, Integer.MAX_VALUE);
	
    /**
     * New pane for configuring the GATKBAMtoPileup node.
     */
    protected GATKBAMtoPileupNodeDialog() {

    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
//    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
//    	createNewGroup("GATK Infile");
//    	DialogComponentFileChooser infile= new DialogComponentFileChooser(INFILE, "infile_variant_filter", JFileChooser.OPEN_DIALOG, false, ".vcf");
////    	gatkf.setBorderTitle("Choose File for filtering");
//    	addDialogComponent(infile);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
//    	gatkf.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
//    	createNewGroup("PED File");
//    	DialogComponentFileChooser ped_file= new DialogComponentFileChooser(BED_FILE, "ped_file", JFileChooser.OPEN_DIALOG, false, ".bed");
////    	gatkf.setBorderTitle("Choose the PED file)");
//    	addDialogComponent(ped_file);
    	
        createNewGroup("Number of Data threads");
        addDialogComponent(new DialogComponentNumber(GATK_NT, "Threads", 1));
        
        createNewGroup("Number of Core threads (for each data thread)");
        addDialogComponent(new DialogComponentNumber(GATK_NCT, "Threads", 1));
    }
}

