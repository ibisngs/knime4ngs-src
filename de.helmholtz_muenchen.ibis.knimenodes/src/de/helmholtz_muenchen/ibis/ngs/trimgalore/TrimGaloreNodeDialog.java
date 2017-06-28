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
package de.helmholtz_muenchen.ibis.ngs.trimgalore;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "TrimGalore" Node.
 * Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
 *
 * 
 * @author Paul Hager
 */
public class TrimGaloreNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring TrimGalore node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
	private SettingsModelBoolean fastqc_enabled;
	private SettingsModelString fastqc_additional_options;
	private SettingsModelIntegerBounded threads;
	private SettingsModelString outfolder_fastqc;
	
	
	

    protected TrimGaloreNodeDialog() {}


	@Override
	public void addToolDialogComponents() {
		
		final int defaultFastQCThreads = TrimGaloreNodeModel.DEFAULT_FASTQC_THREADS;
		final int defaultQuality = TrimGaloreNodeModel.DEFAULT_QUALITY;
		final int defaultStringency = TrimGaloreNodeModel.DEFAULT_STRINGENCY;
		final double defaultErrorRate = TrimGaloreNodeModel.DEFAULT_ERROR_RATE;
		final int defaultLength = TrimGaloreNodeModel.DEFAULT_LENGTH;
		
		final SettingsModelString cutadapt 					= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_CUTADAPT, "");
		final SettingsModelString fastqc 					= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_FASTQC,"");
		final SettingsModelString trimg						= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_TRIMG, "");
		
		fastqc_enabled										= new SettingsModelBoolean(TrimGaloreNodeModel.CFGKEY_FASTQC_ENABLE, false);
		threads 											= new SettingsModelIntegerBounded(TrimGaloreNodeModel.CFGKEY_THREADS_FASTQC, defaultFastQCThreads, 1, Integer.MAX_VALUE);	
		outfolder_fastqc 									= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_OUTFOLDER_FASTQC,"");    	
		
		final SettingsModelString outfolder_trimg 			= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_OUTFOLDER_TRIMG,"");
		
		final SettingsModelString adapter 					= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_ADAPTER,"");
		final SettingsModelString adapter2 					= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_ADAPTER2,"");
		final SettingsModelString preset_adapter 			= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_PRESET_ADAPTER,"");
		
		final SettingsModelIntegerBounded quality 			= new SettingsModelIntegerBounded(TrimGaloreNodeModel.CFGKEY_QUALITY, defaultQuality, 0, Integer.MAX_VALUE);
		final SettingsModelIntegerBounded stringency 		= new SettingsModelIntegerBounded(TrimGaloreNodeModel.CFGKEY_STRINGENCY, defaultStringency, 1, Integer.MAX_VALUE);
		final SettingsModelDoubleBounded error_rate 		= new SettingsModelDoubleBounded(TrimGaloreNodeModel.CFGKEY_ERROR_RATE, defaultErrorRate, 0, 1);
		
		final SettingsModelBoolean gzip 					= new SettingsModelBoolean(TrimGaloreNodeModel.CFGKEY_GZIP, true);
		
		//final SettingsModelIntegerBounded max_length 		= new SettingsModelIntegerBounded(TrimGaloreNodeModel.CFGKEY_MAX_LENGTH, 0, 0, Integer.MAX_VALUE);
		final SettingsModelIntegerBounded length 			= new SettingsModelIntegerBounded(TrimGaloreNodeModel.CFGKEY_LENGTH, defaultLength, 0, Integer.MAX_VALUE);
		
		final SettingsModelString additional_options 		= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_ADDITIONAL_OPTIONS, "");
		fastqc_additional_options							= new SettingsModelString(TrimGaloreNodeModel.CFGKEY_FASTQC_ADDITIONAL_OPTIONS, "");

		
		boolean enableFastQC = fastqc_enabled.getBooleanValue();
		
		// Preferences page settings
		addPrefPageSetting(fastqc, IBISKNIMENodesPlugin.FASTQC);
    	addPrefPageSetting(trimg, IBISKNIMENodesPlugin.TRIMG);
    	addPrefPageSetting(cutadapt, IBISKNIMENodesPlugin.CUTADAPT);
		
    	
		// Primary TrimGalore options tab
    	setDefaultTabTitle("Trim Galore Options");
    	
    	createNewGroup("Folder for trim galore output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder_trimg, "his_id_TrimGalore_OUT_TRIMG", JFileChooser.OPEN_DIALOG, true));
    	
    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentBoolean(gzip, "Gzip ouput"));
    	
    	addDialogComponent(new DialogComponentString(adapter, "Adapter sequence"));
    	addDialogComponent(new DialogComponentString(adapter2, "Adapter sequence for read file 2 (paired-end)"));
    	addDialogComponent(new DialogComponentStringSelection(preset_adapter, "Preset adapter", "AUTOMATIC", "illumina", "nextera", "small_rna"));
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(quality, "Quality", 1));
    	addDialogComponent(new DialogComponentNumber(stringency, "Stringency", 1));
    	addDialogComponent(new DialogComponentNumber(error_rate, "Error rate", 0.1));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	//addDialogComponent(new DialogComponentNumberEdit(max_length, "Max length"));
    	addDialogComponent(new DialogComponentNumber(length, "Length", 1));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentString(additional_options, "Additional Trim Galore options"));
    	
    	// Enable/disable preset_adapter based on if custom adapter has been defined or not
    	adapter.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				if(!adapter.getStringValue().equals("")){
					preset_adapter.setEnabled(false);
				} else {
					preset_adapter.setEnabled(true);
				}
			}
		});
    	
    	
    	
    	// FastQC Options
    	createNewTab("FastQC Options");
    	addDialogComponent(new DialogComponentBoolean(fastqc_enabled, "Execute FastQC after Trim Galore"));
    	
    	outfolder_fastqc.setEnabled(enableFastQC);
    	threads.setEnabled(enableFastQC);
    	fastqc_additional_options.setEnabled(enableFastQC);
    	
    	createNewGroup("Folder for fastqc output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder_fastqc, "his_id_TrimGalore_OUT_FASTQC", JFileChooser.OPEN_DIALOG, true));

    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1));
    	addDialogComponent(new DialogComponentString(fastqc_additional_options, "Additional FastQC options"));
    	
    	// Makes the fastqc options tab dependant on the fastqc checkbox
		fastqc_enabled.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(final ChangeEvent e) {
				
				boolean enableFastQC = fastqc_enabled.getBooleanValue();
				
				outfolder_fastqc.setEnabled(enableFastQC);
				threads.setEnabled(enableFastQC);
				fastqc_additional_options.setEnabled(enableFastQC);
			}
		});
	}
	
	public void onOpen() {
		super.onOpen();
		
		boolean enableFastQC = fastqc_enabled.getBooleanValue();
		
		outfolder_fastqc.setEnabled(enableFastQC);
    	threads.setEnabled(enableFastQC);
    	fastqc_additional_options.setEnabled(enableFastQC);
	}
}

