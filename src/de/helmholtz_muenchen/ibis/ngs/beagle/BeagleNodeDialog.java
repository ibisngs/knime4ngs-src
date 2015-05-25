package de.helmholtz_muenchen.ibis.ngs.beagle;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "Beagle" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tanzeem
 */
public class BeagleNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring Beagle node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	private final SettingsModelString BEAGLE = new SettingsModelString(BeagleNodeModel.CFGKEY_BEAGLE_PATH, "");
    private final SettingsModelString INFILE = new SettingsModelString(BeagleNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString REF_VCF = new SettingsModelString(BeagleNodeModel.CFGKEY_REF_VCF, "");
    private final SettingsModelString PED_FILE = new SettingsModelString(BeagleNodeModel.CFGKEY_PED_FILE, "");
	final SettingsModelString param = new SettingsModelString(BeagleNodeModel.CFGKEY_PARAM, "");	//tool selection

    
    protected BeagleNodeDialog() {
        super();
        
        createNewGroup("Path to BEAGLE jar file");
    	DialogComponentFileChooser beagle_jar= new DialogComponentFileChooser(BEAGLE, "beagle", 0, ".jar");
//    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(beagle_jar);
    	
    	createNewGroup("BEAGLE Input");
    	DialogComponentFileChooser infile= new DialogComponentFileChooser(INFILE, "infile_variant_filter", 0, ".vcf");
    	infile.setBorderTitle("Select input (disabled if file available from previous node)");
    	addDialogComponent(infile);
    	
    	createNewGroup("Reference phased file");
    	DialogComponentFileChooser ref_vcf= new DialogComponentFileChooser(REF_VCF, "ref_genome_phased", 0, ".vcf");
    	ref_vcf.setBorderTitle("Select reference (if available)");
    	addDialogComponent(ref_vcf);
    	
    	createNewGroup("PED File");
    	DialogComponentFileChooser ped_file= new DialogComponentFileChooser(PED_FILE, "ped_file", /*JFileChooser.OPEN_DIALOG, false,*/0, ".ped");
//    	gatkf.setBorderTitle("Choose the PED file)");
    	addDialogComponent(ped_file);
    	 	
    	//tool selection
        createNewGroup("Input parameter selection");
        addDialogComponent(new DialogComponentStringSelection(param, "Parameter", BeagleNodeModel.GENOTYPE_PARAMS));
        
                    
    }
}

