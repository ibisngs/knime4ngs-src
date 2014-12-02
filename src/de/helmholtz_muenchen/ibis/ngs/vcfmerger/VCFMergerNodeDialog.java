package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VCFMerger" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeDialog extends DefaultNodeSettingsPane {

	  private final SettingsModelString INFOLDER 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_INFOLDER, "");
	  private final SettingsModelString REGEX 					= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REGEX, "");
	  private final SettingsModelString OUTFOLDER 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");
	  private final SettingsModelString GATK 					= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GATK, "");
	  private final SettingsModelString REF_GENOME 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REF_GENOME, "");
	  private final SettingsModelString GENOTYPEMERGEOPTION	= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
	  private final SettingsModelString OUTFILETAG			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFILETAG, "");
	
    /**
     * New pane for configuring the VCFMerger node.
     */
    protected VCFMergerNodeDialog() {
    		
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
//    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
 	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
//    	gatkf.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
    	
    	createNewGroup("Folder in which the search is performed");
    	addDialogComponent(new DialogComponentFileChooser(INFOLDER, "Folder in which the search is performed", 0,true,""));
    	createNewGroup("Outfolder");
    	addDialogComponent(new DialogComponentFileChooser(OUTFOLDER, "Outfolder", 0,true,""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(REGEX, "File Suffix of VCF Files to merge"));
    	addDialogComponent(new DialogComponentString(OUTFILETAG, "Outfile Tag of merged VCF File"));
    	addDialogComponent(new DialogComponentStringSelection(GENOTYPEMERGEOPTION, "Genotype Merge Strategy","UNSORTED","UNIQUIFY","REQUIRE_UNIQUE"));
    	
    }
}

