package de.helmholtz_muenchen.ibis.ngs.vep;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
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
	
	//VEP options tab
	private final SettingsModelString veppl = new SettingsModelString(VEPNodeModel.CFGKEY_VEP_PL,"-");
	private final SettingsModelString fasta = new SettingsModelString(VEPNodeModel.CFGKEY_FASTA_FILE,VEPNodeModel.DEF_CACHE_DIR);
	private final SettingsModelString outfolder = new SettingsModelString(VEPNodeModel.CFGKEY_OUTFOLDER,"");
	private final SettingsModelBoolean overwrite = new SettingsModelBoolean(VEPNodeModel.CFGKEY_OVERWRITE,false);
	private final SettingsModelBoolean use_cache = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_CACHE,true);
	private final SettingsModelInteger forks = new SettingsModelInteger(VEPNodeModel.CFGKEY_FORKS,1);
	private final SettingsModelBoolean coding_only = new SettingsModelBoolean(VEPNodeModel.CFGKEY_CODING_ONLY, true);
	private final SettingsModelString transcript_set = new SettingsModelString(VEPNodeModel.CFGKEY_TRANSCRIPT_SET,"GENCODE Basic");

	//advanced tab
	private final SettingsModelString stats_type = new SettingsModelString(VEPNodeModel.CFGKEY_STATS_TYPE,"html");
	private final SettingsModelString cache_dir = new SettingsModelString(VEPNodeModel.CFGKEY_CACHE_DIR, VEPNodeModel.DEF_CACHE_DIR);
	private final SettingsModelString plugin_dir = new SettingsModelString(VEPNodeModel.CFGKEY_PLUGIN_DIR, VEPNodeModel.DEF_PLUGIN_DIR);
	private final SettingsModelString further_plugins = new SettingsModelString(VEPNodeModel.CFGKEY_FURTHER_PLUGINS, "");
	
	//LOFTEE tab
	private final SettingsModelBoolean use_loftee = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_LOFTEE,false);
	private final SettingsModelString human_ancestor = new SettingsModelString(VEPNodeModel.CFGKEY_HUMAN_ANCESTOR,"-");
	private final SettingsModelString conservation_file = new SettingsModelString(VEPNodeModel.CFGKEY_CONSERVATION_FILE,"-");
	private final SettingsModelString samtools_path = new SettingsModelString(VEPNodeModel.CFGKEY_SAMTOOLS_PATH, "-");
	
    protected VEPNodeDialog() {

    	addPrefPageSetting(veppl, IBISKNIMENodesPlugin.VEP);
    	addPrefPageSetting(samtools_path, IBISKNIMENodesPlugin.SAMTOOLS);
    	
    	createNewGroup("Path to variant_effect_predictor.pl");
    	addDialogComponent(new DialogComponentFileChooser(veppl, "his_id_VEP_VEPPL", 0, ".pl"));
    	
    	createNewGroup("Path to fasta file");
    	addDialogComponent(new DialogComponentFileChooser(fasta, "his_id_FASTA", 0, ".fa|.fasta"));
    	
    	createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_VEP_OUT", 0, true));
    	addDialogComponent(new DialogComponentBoolean(overwrite, "Overwrite, if output files exist?"));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(use_cache, "Use cache?"));
    	addDialogComponent(new DialogComponentNumber(forks, "Number of forks",1));
    	addDialogComponent(new DialogComponentBoolean(coding_only, "Coding only?"));
    	addDialogComponent(new DialogComponentStringSelection(transcript_set,"Choose transcript set", VEPNodeModel.TRANSCRIPT_SETS));

    	createNewTab("Advanced");
    	addDialogComponent(new DialogComponentStringSelection(stats_type,"Choose type of statistic file", VEPNodeModel.STAT_TYPES));
    	
    	createNewGroup("Cache directory");
    	addDialogComponent(new DialogComponentFileChooser(cache_dir, "his_id_VEP_CACHEDIR", 0, true));
    	
    	createNewGroup("Plugin directory");
    	addDialogComponent(new DialogComponentFileChooser(plugin_dir, "his_id_VEP_PLUGINDIR", 0, true));
    	
    	createNewGroup("Specify further plugins/options");
    	addDialogComponent(new DialogComponentString(further_plugins, "Plugins/options"));
    	
    	createNewTab("LOFTEE");
    	addDialogComponent(new DialogComponentBoolean(use_loftee, "Use LOFTEE?"));
    	
    	createNewGroup("Path to human_ancestor.fa.gz");
    	addDialogComponent(new DialogComponentFileChooser(human_ancestor, "his_id_VEP_HUMANANCESTOR", 0, ".fa.gz"));
    	
    	createNewGroup("Path to phylocsf.sql");
    	addDialogComponent(new DialogComponentFileChooser(conservation_file, "his_id_VEP_CONSERVATIONFILE", 0, ".sql"));
    	
    	createNewGroup("Samtools PATH");
    	addDialogComponent(new DialogComponentFileChooser(samtools_path, "his_id_VEP_SAMTOOLSPATH", 0, ""));
  
//    	use_loftee.addChangeListener(new ChangeListener () {
//
//			@Override
//			public void stateChanged(ChangeEvent e) {
//				updatePrefs();
//				
//			}
//    		
//    	});
    }

//	@Override
//	protected void updatePrefs() {
//		if(usePrefPage.getBooleanValue()) {
//	    	String vep_path = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("variant_effect_predictor.pl");
//	    	if(vep_path != null && !vep_path.equals("")) {
//	    		veppl.setStringValue(vep_path);
//	    		veppl.setEnabled(false);
//			} else {
//				veppl.setEnabled(true);
//			}
//	    	if(use_loftee.getBooleanValue()){
//		    	String toolPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("samtools");
//		    	if(toolPath != null && !toolPath.equals("")) {
//		    		samtools_path.setStringValue(toolPath);
//		    		samtools_path.setEnabled(false);
//		    	} else {
//		    		samtools_path.setEnabled(true);
//		    	}
//	    	}
//		} else {
//			veppl.setEnabled(true);
//			samtools_path.setEnabled(true);
//		}
//	}
}

