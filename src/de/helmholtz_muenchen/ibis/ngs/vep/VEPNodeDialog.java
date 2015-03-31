package de.helmholtz_muenchen.ibis.ngs.vep;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "VEP" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class VEPNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VEP node.
     */
	
	//VEP tab
	private final SettingsModelString vcfin = new SettingsModelString(VEPNodeModel.CFGKEY_VCF_INFILE,"-");
	private final SettingsModelString veppl = new SettingsModelString(VEPNodeModel.CFGKEY_VEP_PL,"-");
	private final SettingsModelString stats_type = new SettingsModelString(VEPNodeModel.CFGKEY_STATS_TYPE,"html");
	private final SettingsModelBoolean coding_only = new SettingsModelBoolean(VEPNodeModel.CFGKEY_CODING_ONLY, true);
	private final SettingsModelInteger forks = new SettingsModelInteger(VEPNodeModel.CFGKEY_FORKS,1);
	private final SettingsModelBoolean use_cache = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_CACHE,true);
	private final SettingsModelString cache_dir = new SettingsModelString(VEPNodeModel.CFGKEY_CACHE_DIR, VEPNodeModel.DEF_CACHE_DIR);
	private final SettingsModelString plugin_dir = new SettingsModelString(VEPNodeModel.CFGKEY_PLUGIN_DIR, VEPNodeModel.DEF_CACHE_DIR);

	//LOFTEE tab
	private final SettingsModelBoolean use_loftee = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_LOFTEE,false);
	private final SettingsModelString human_ancestor = new SettingsModelString(VEPNodeModel.CFGKEY_HUMAN_ANCESTOR,"-");
	private final SettingsModelString conservation_file = new SettingsModelString(VEPNodeModel.CFGKEY_CONSERVATION_FILE,"-");
	private final SettingsModelString samtools_path = new SettingsModelString(VEPNodeModel.CFGKEY_SAMTOOLS_PATH, "-");

	
    protected VEPNodeDialog() {

    	super();
    	
    	createNewGroup("Path to VCF file");
    	addDialogComponent(new DialogComponentFileChooser(vcfin, "his_id_VEP_VCFIN", 0, ".vcf"));
    	
    	createNewGroup("Path to variant_effect_predictor.pl");
    	addDialogComponent(new DialogComponentFileChooser(veppl, "his_id_VEP_VEPPL", 0, ".pl"));
    	
    	createNewGroup("Cache usage");
    	addDialogComponent(new DialogComponentBoolean(use_cache, "Use cache?"));
    	addDialogComponent(new DialogComponentFileChooser(cache_dir, "his_id_VEP_CACHEDIR", 0, true));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(coding_only, "Coding only?"));
    	addDialogComponent(new DialogComponentButtonGroup(stats_type, false, "Choose type of statistic file", VEPNodeModel.STAT_TYPES));
    	addDialogComponent(new DialogComponentNumber(forks, "Number of forks",1));
    	
    	createNewGroup("Plugin directory");
    	addDialogComponent(new DialogComponentFileChooser(plugin_dir, "his_id_VEP_PLUGINDIR", 0, true));

    	createNewTab("LOFTEE");
    	addDialogComponent(new DialogComponentBoolean(use_loftee, "Use LOFTEE?"));
    	
    	createNewGroup("Path to human_ancestor.fa.gz");
    	addDialogComponent(new DialogComponentFileChooser(human_ancestor, "his_id_VEP_HUMANANCESTOR", 0, ".fa.gz"));
    	
    	createNewGroup("Path to phylocsf.sql");
    	addDialogComponent(new DialogComponentFileChooser(conservation_file, "his_id_VEP_CONSERVATIONFILE", 0, ".sql"));
    	
    	createNewGroup("Samtools PATH");
    	addDialogComponent(new DialogComponentFileChooser(samtools_path, "his_id_VEP_SAMTOOLSPATH", 0, true));
    }
}

