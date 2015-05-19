package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKPhaseByTransmission" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKPhaseByTransmissionNodeDialog extends DefaultNodeSettingsPane {

	
	    private final SettingsModelString GATK = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_GATK_PATH, "");
	    private final SettingsModelString INFILE = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_INFILE, "");
	    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_REF_GENOME, "");
	    private final SettingsModelString PED_FILE = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_PED_FILE, "");
	    private final SettingsModelIntegerBounded GATK_JAVA_MEMORY = new SettingsModelIntegerBounded(GATKPhaseByTransmissionNodeModel.CFGKEY_JAVAMEMORY, GATKPhaseByTransmissionNodeModel.DEF_NUM_JAVAMEMORY, GATKPhaseByTransmissionNodeModel.MIN_NUM_JAVAMEMORY, GATKPhaseByTransmissionNodeModel.MAX_NUM_JAVAMEMORY);
	    private final SettingsModelString DENOVO_PRIOR = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_DENOVOPRIOR, "1.0E-8");
	
	
    /**
     * New pane for configuring the GATKPhaseByTransmission node.
     */
    protected GATKPhaseByTransmissionNodeDialog() {

    	
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", 0, ".jar");
//    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
    	createNewGroup("GATK Infile");
    	DialogComponentFileChooser infile= new DialogComponentFileChooser(INFILE, "infile_variant_filter", 0, ".vcf");
    	infile.setBorderTitle("Select input (disabled if file available from previous node)");
    	addDialogComponent(infile);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", 0, ".fa", ".fasta", ".FASTA");
    	ref_genome.setBorderTitle("Select reference (disabled if file available as a flow variable)");
    	addDialogComponent(ref_genome);
    	
    	createNewGroup("PED File");
    	DialogComponentFileChooser ped_file= new DialogComponentFileChooser(PED_FILE, "ped_file", /*JFileChooser.OPEN_DIALOG, false,*/0, ".ped");
//    	gatkf.setBorderTitle("Choose the PED file)");
    	addDialogComponent(ped_file);
    	
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(DENOVO_PRIOR, "De Novo Prior"));
    	
        createNewGroup("Java Memory");
        addDialogComponent(new DialogComponentNumber(GATK_JAVA_MEMORY, "Java Memory (GB)", 1));
    	
    }
}

