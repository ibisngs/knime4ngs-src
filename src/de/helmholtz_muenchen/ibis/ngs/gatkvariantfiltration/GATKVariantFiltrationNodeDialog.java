package de.helmholtz_muenchen.ibis.ngs.gatkvariantfiltration;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKVariantFiltration" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKVariantFiltrationNodeDialog extends DefaultNodeSettingsPane {

	
    private final SettingsModelString GATK = new SettingsModelString(GATKVariantFiltrationNodeModel.CFGKEY_GATK_PATH, "");
//    private final SettingsModelString INFILE = new SettingsModelString(GATKVariantFiltrationNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKVariantFiltrationNodeModel.CFGKEY_REF_GENOME, "");
    
	private final SettingsModelOptionalString QUAL= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_QUAL, "< 50.0",true);
	private final SettingsModelOptionalString DP= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_DP, "< 20.0",true);
	private final SettingsModelOptionalString AD= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_AD, "< 5.0",false);
	private final SettingsModelOptionalString QD= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_QD, "< 2.0",false);
	private final SettingsModelOptionalString FS= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FS, "> 60.0",false);
	private final SettingsModelOptionalString MQ= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_MQ, "< 40.0",false);
	private final SettingsModelOptionalString HS= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_HS, "> 13.0",false);
	private final SettingsModelOptionalString MQR= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_MQR, "< -12.5",false);
	private final SettingsModelOptionalString RPR= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_RPR, "< -8.0",false);
	private final SettingsModelOptionalString FilterString = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FilterString, "-",false);
	private final SettingsModelOptionalString FilterName = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FilterName, "-",false);
	private final SettingsModelIntegerBounded memory_usage = new SettingsModelIntegerBounded(GATKVariantFiltrationNodeModel.CFGKEY_JAVAMEMORY, GATKVariantFiltrationNodeModel.DEF_NUM_JAVAMEMORY, GATKVariantFiltrationNodeModel.MIN_NUM_JAVAMEMORY, GATKVariantFiltrationNodeModel.MAX_NUM_JAVAMEMORY);


	
    /**
     * New pane for configuring the GATKVariantFiltration node.
     */
    protected GATKVariantFiltrationNodeDialog() {

    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
//    	createNewGroup("GATK Infile");
//    	DialogComponentFileChooser infile= new DialogComponentFileChooser(INFILE, "infile_variant_filter", JFileChooser.OPEN_DIALOG, false, ".vcf");
//    	gatkf.setBorderTitle("Choose File for filtering");
//    	addDialogComponent(infile);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	gatkf.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
        //#Memory
        createNewGroup("");
        addDialogComponent(new DialogComponentNumber(memory_usage, "Java Memory (GB)", 1));
    	
    	createNewTab("Filter Options");
    	
    	createNewGroup("Filter Name");
    	addDialogComponent(new DialogComponentOptionalString(FilterName, "Filter Name"));
    	
    	createNewGroup("Discard all variants that match any of the following conditions:");
    	addDialogComponent(new DialogComponentOptionalString(QUAL, "Quality Cutoff (QUAL)"));
    	addDialogComponent(new DialogComponentOptionalString(DP, "Coverage Cutoff (DP)"));
    	addDialogComponent(new DialogComponentOptionalString(AD, "Allelic depth"));
    	addDialogComponent(new DialogComponentOptionalString(QD, "Variant Confidence/Quality by Depth (QD)"));
    	addDialogComponent(new DialogComponentOptionalString(FS, "Strand Bias (FS)"));
    	addDialogComponent(new DialogComponentOptionalString(MQ, "RMS Mapping Quality (MQ)"));
    	addDialogComponent(new DialogComponentOptionalString(HS, "HaplotypeScore"));
    	addDialogComponent(new DialogComponentOptionalString(MQR, "MappingQualityRankSum"));
    	addDialogComponent(new DialogComponentOptionalString(RPR, "ReadPosRankSum"));
    	addDialogComponent(new DialogComponentOptionalString(FilterString, "Additional Filter Options",50));

    }
}

