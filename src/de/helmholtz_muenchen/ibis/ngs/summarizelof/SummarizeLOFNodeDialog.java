package de.helmholtz_muenchen.ibis.ngs.summarizelof;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "SummarizeLOF" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class SummarizeLOFNodeDialog extends DefaultNodeSettingsPane {

	
	private final SettingsModelString pedfile = new SettingsModelString(SummarizeLOFNodeModel.CFGKEY_PED_FILE,"-");
	private final SettingsModelString vcfin = new SettingsModelString(SummarizeLOFNodeModel.CFGKEY_VCF_INFILE,"-");
	
	
    /**
     * New pane for configuring the SummarizeLOF node.
     */
    protected SummarizeLOFNodeDialog() {
    	
    	createNewGroup("VCF File");
    	addDialogComponent(new DialogComponentFileChooser(vcfin, "his_id_SummarizeLOF_VCFIN", 0, ".vcf"));
    	createNewGroup("PED File");
    	addDialogComponent(new DialogComponentFileChooser(pedfile, "his_id_SummarizeLOF_PED", 0, ".ped")); 	
    }
}

