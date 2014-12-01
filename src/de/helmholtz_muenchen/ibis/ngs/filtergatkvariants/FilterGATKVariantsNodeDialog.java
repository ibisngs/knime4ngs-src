package de.helmholtz_muenchen.ibis.ngs.filtergatkvariants;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "FilterGATKVariants" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class FilterGATKVariantsNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString m_PINDEL_INFILE = new SettingsModelString(FilterGATKVariantsNodeModel.CFGKEY_PINDEL_INFILE, "");
	private final SettingsModelString m_GATK_INFILE = new SettingsModelString(FilterGATKVariantsNodeModel.CFGKEY_GATK_INFILE, "");
	
	
    /**
     * New pane for configuring the FilterGATKVariants node.
     */
    protected FilterGATKVariantsNodeDialog() {

    	createNewGroup("Pindel Infile");
    	addDialogComponent(new DialogComponentFileChooser(m_PINDEL_INFILE, "pindel_infile", 0, ".vcf"));
    	createNewGroup("GATK Infile");
    	addDialogComponent(new DialogComponentFileChooser(m_GATK_INFILE, "gatk_pindel_infile", 0, ".vcf"));
    	
    }
}

