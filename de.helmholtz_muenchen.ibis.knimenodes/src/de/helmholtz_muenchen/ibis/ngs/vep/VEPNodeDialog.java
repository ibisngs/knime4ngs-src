/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.vep;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;


/**
 * <code>NodeDialog</code> for the "VEP" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class VEPNodeDialog extends HTExecutorNodeDialog {
	
	public VEPNodeDialog(){
		super(VCFCell.TYPE.getPreferredValueClass(), 0);
	}

	private SettingsModelBoolean my_overwrite;
    
    public void addToolDialogComponents() {

    	//VEP options tab
    	final SettingsModelString veppl = new SettingsModelString(VEPNodeModel.CFGKEY_VEP_PL,"-");
    	final SettingsModelString fasta = new SettingsModelString(VEPNodeModel.CFGKEY_FASTA_FILE,VEPNodeModel.DEF_CACHE_DIR);
        final SettingsModelString m_out_format = new SettingsModelString(VEPNodeModel.CFGKEY_OUT_FORMAT,VEPNodeModel.DEF_FORMAT);
    	final SettingsModelString outfolder = new SettingsModelString(VEPNodeModel.CFGKEY_OUTFOLDER,"");
    	final SettingsModelBoolean overwrite = new SettingsModelBoolean(VEPNodeModel.CFGKEY_OVERWRITE,false);
    	my_overwrite = overwrite;
    	final SettingsModelBoolean use_cache = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_CACHE,true);
    	final SettingsModelOptionalString further_options = new SettingsModelOptionalString(VEPNodeModel.CFGKEY_FURTHER_OPTIONS, "", false);
    	final SettingsModelInteger forks = new SettingsModelInteger(VEPNodeModel.CFGKEY_FORKS,1);
    	final SettingsModelInteger buffer = new SettingsModelInteger(VEPNodeModel.CFGKEY_BUFFER, 5000);
    	final SettingsModelBoolean coding_only = new SettingsModelBoolean(VEPNodeModel.CFGKEY_CODING_ONLY, true);
    	final SettingsModelString transcript_set = new SettingsModelString(VEPNodeModel.CFGKEY_TRANSCRIPT_SET,"GENCODE Basic");
        final SettingsModelBoolean sift = new SettingsModelBoolean(VEPNodeModel.CFGKEY_SIFT, true);
        final SettingsModelBoolean polyphen = new SettingsModelBoolean(VEPNodeModel.CFGKEY_POLYPHEN, true);
        final SettingsModelBoolean symbol = new SettingsModelBoolean(VEPNodeModel.CFGKEY_SYMBOL, true);
        final SettingsModelBoolean biotype = new SettingsModelBoolean(VEPNodeModel.CFGKEY_BIOTYPE, true);
        
    	//advanced tab
    	final SettingsModelString stats_type = new SettingsModelString(VEPNodeModel.CFGKEY_STATS_TYPE,"html");
    	final SettingsModelString cache_dir = new SettingsModelString(VEPNodeModel.CFGKEY_CACHE_DIR, VEPNodeModel.DEF_CACHE_DIR);
    	final SettingsModelString plugin_dir = new SettingsModelString(VEPNodeModel.CFGKEY_PLUGIN_DIR, VEPNodeModel.DEF_PLUGIN_DIR);
    	final SettingsModelOptionalString further_plugins = new SettingsModelOptionalString(VEPNodeModel.CFGKEY_FURTHER_PLUGINS, "", false);
    	
    	//LOFTEE tab
    	final SettingsModelBoolean use_loftee = new SettingsModelBoolean(VEPNodeModel.CFGKEY_USE_LOFTEE,false);
    	final SettingsModelString human_ancestor = new SettingsModelString(VEPNodeModel.CFGKEY_HUMAN_ANCESTOR,"-");
    	final SettingsModelString conservation_file = new SettingsModelString(VEPNodeModel.CFGKEY_CONSERVATION_FILE,"-");
    	final SettingsModelString samtools_path = new SettingsModelString(VEPNodeModel.CFGKEY_SAMTOOLS_PATH, "-");
    	
    	addPrefPageSetting(veppl, IBISKNIMENodesPlugin.VEP);
    	addPrefPageSetting(samtools_path, IBISKNIMENodesPlugin.SAMTOOLS);
    	
    	setDefaultTabTitle("Annotation");
    	addDialogComponent(new DialogComponentBoolean(coding_only, "Coding only?"));
    	addDialogComponent(new DialogComponentStringSelection(transcript_set,"Choose transcript set", VEPNodeModel.TRANSCRIPT_SETS));
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(sift, "SIFT"));
    	addDialogComponent(new DialogComponentBoolean(polyphen, "PolyPhen"));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentBoolean(symbol, "Add gene symbol?"));
    	addDialogComponent(new DialogComponentBoolean(biotype, "Add biotype of transcript or regulatory feature?"));
    	addDialogComponent(new DialogComponentOptionalString(further_options, "Further flags"));
    	
    	createNewTab("Performance");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(forks, "Number of forks",1));
    	addDialogComponent(new DialogComponentNumber(buffer, "Buffer size",1000));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentBoolean(use_cache, "Use cache?"));
    	createNewGroup("Cache directory");
    	addDialogComponent(new DialogComponentFileChooser(cache_dir, "his_id_VEP_CACHEDIR", 0, true));
    	
    	createNewGroup("Path to fasta file");
    	addDialogComponent(new DialogComponentFileChooser(fasta, "his_id_FASTA", 0, ".fa|.fasta"));

    	
    	createNewTab("Output");
    	addDialogComponent(new DialogComponentStringSelection(m_out_format, "Choose format of annotation file" ,VEPNodeModel.OUT_FORMATS));
    	addDialogComponent(new DialogComponentStringSelection(stats_type,"Choose format of statistics file", VEPNodeModel.STAT_TYPES));
    	
    	createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_VEP_OUT", 0, true));
    	addDialogComponent(new DialogComponentBoolean(overwrite, "Overwrite, if output files exist?"));
    	
    	createNewTab("Plugins");
    	
    	createNewGroup("Plugin directory");
    	addDialogComponent(new DialogComponentFileChooser(plugin_dir, "his_id_VEP_PLUGINDIR", 0, true));
    	
    	addDialogComponent(new DialogComponentOptionalString(further_plugins, "Further plugins"));
    	
    	createNewTab("LOFTEE");
    	addDialogComponent(new DialogComponentBoolean(use_loftee, "Use LOFTEE?"));
    	
    	createNewGroup("Path to human_ancestor.fa.gz");
    	addDialogComponent(new DialogComponentFileChooser(human_ancestor, "his_id_VEP_HUMANANCESTOR", 0, ".fa.gz"));
    	
    	createNewGroup("Path to phylocsf.sql");
    	addDialogComponent(new DialogComponentFileChooser(conservation_file, "his_id_VEP_CONSERVATIONFILE", 0, ".sql"));
    	
    	use_cache.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(use_cache.getBooleanValue()) {
					cache_dir.setEnabled(true);
				} else {
					cache_dir.setEnabled(false);
				}
			}
    	});
    }
    
    
    
    public void onOpen() {
    	super.onOpen();
    	boolean use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
    	
    	if(use_hte) {
    		my_overwrite.setBooleanValue(true);
    		my_overwrite.setEnabled(false);
    	} else {
//    		my_overwrite.setBooleanValue(false);
    		my_overwrite.setEnabled(true);
    	}
    }

}

