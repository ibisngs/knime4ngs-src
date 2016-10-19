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
package de.helmholtz_muenchen.ibis.ngs.bcftools;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Bcftools" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class BcftoolsNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the Bcftools node.
     */
    protected BcftoolsNodeDialog() {
    }
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString bcfpath = new SettingsModelString(BcftoolsNodeModel.CFGKEY_PATH2BCFTOOLS, "");
    	final SettingsModelString bcfmethod = new SettingsModelString(BcftoolsNodeModel.CFGKEY_BCFMETHOD,"");
    	final SettingsModelString vcfsampleheader = new SettingsModelString(BcftoolsNodeModel.CFGKEY_VCFSAMPLEHEADER, "");
    	vcfsampleheader.setEnabled(false);
    	
    	addPrefPageSetting(bcfpath, IBISKNIMENodesPlugin.BCFTOOLS);
    	
    	//Concat
    	final SettingsModelBoolean concat_overlap 	= new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_CONCAT_OVERLAPS, false);
    	final SettingsModelString concat_outfile_type = new SettingsModelString(
    			BcftoolsNodeModel.CFGKEY_CONCAT_OUTFILE_TYPE, "uncompressed VCF");
    	
    	//Extra
    	final SettingsModelOptionalString furtherOptions = new SettingsModelOptionalString(
    			BcftoolsNodeModel.CFGKEY_FURTHER_OPTIONS,"",false);

    	
//    	//Call
//    	final SettingsModelBoolean ifbedfile =  new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_IFBEDFILE, false);
//    	final SettingsModelBoolean ifsamplelist =  new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_IFSAMPLELIST, false);
//    	final SettingsModelString bedfile = new SettingsModelString(BcftoolsNodeModel.CFGKEY_BEDFILE, "");
//    	final SettingsModelString samplelist = new SettingsModelString(BcftoolsNodeModel.CFGKEY_SAMPLELIST, "");
//    	bedfile.setEnabled(false);
//    	samplelist.setEnabled(false);
//    	final SettingsModelBoolean snpcalling = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_SNPCALLING, false);	
//    	final SettingsModelBoolean outvariantsonly = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_OUTVARIANTSONLY, false);
//    	final SettingsModelString constrainedcalling = new SettingsModelString(BcftoolsNodeModel.CFGKEY_CONSTRAINEDCALLING,"No constrains");
//    	final SettingsModelString call_outfile_type = new SettingsModelString(
//    			BcftoolsNodeModel.CFGKEY_CALL_OUTFILE_TYPE, "");
    	
    	//Main Tab 
//    	createNewGroup("Path to Bcftools");
//    	addDialogComponent(new DialogComponentFileChooser(bcfpath,"his_bcft_ID5",""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentStringSelection(bcfmethod,"Select method","index", "concat","reheader","stats"));
    	createNewGroup("Select file that includes the new sample names");
    	addDialogComponent(new DialogComponentFileChooser(vcfsampleheader,"hisID_bcft_vcfsampleheader",""));
    	addDialogComponent(new DialogComponentOptionalString(furtherOptions,"Further Parameters"));
    	
    	//Concat
    	createNewTab("Concat");
    	createNewGroup("Concat");
    	addDialogComponent(new DialogComponentBoolean(concat_overlap, "First coordinate of the next file can precede last record of the current file"));
    	addDialogComponent(new DialogComponentStringSelection(concat_outfile_type,"Output format:","uncompressed VCF","compressed VCF","uncompressed BCF","compressed BCF"));
    	
//    	createNewTab("Call");
//    	createNewGroup("Call");
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(ifbedfile, "List of sites (chr pos) or regions (BED) to output:"));
//    	addDialogComponent(new DialogComponentFileChooser(bedfile, "his_bcft_ID1",""));
//    	setHorizontalPlacement(false);
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(ifsamplelist, "List of samples to use:"));
//    	addDialogComponent(new DialogComponentFileChooser(samplelist, "his_bcft_ID2",""));
//    	setHorizontalPlacement(false);
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
//    			BcftoolsNodeModel.CFGKEY_KEEPALLELES, false), "Keep all possible alternate alleles at variant sites"));
//    	setHorizontalPlacement(false);
//    	addDialogComponent(new DialogComponentStringSelection(call_outfile_type,"Output format:","uncompressed VCF","compressed VCF","uncompressed BCF","compressed BCF"));
//
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(snpcalling, "Call SNPs"));
//    	setHorizontalPlacement(false);
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(outvariantsonly, "Output potential variant sites only"));
//    	setHorizontalPlacement(false);
//
//
//    	addDialogComponent(new DialogComponentStringSelection(
//    			constrainedcalling,
//    			"Constrained calling","No constrains","pair", "trioauto", "trioxd","trioxs"));


 
    	
    	
    	
    	//Main method chooser
    	bcfmethod.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(bcfmethod.getStringValue().equals("reheader")){
					vcfsampleheader.setEnabled(true);
				}else{
					vcfsampleheader.setEnabled(false);
				}
				
				}
		});
    	
    	
    	

//    	ifbedfile.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				bedfile.setEnabled(ifbedfile.getBooleanValue());
//				}
//		});
//
//    	ifsamplelist.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				samplelist.setEnabled(ifsamplelist.getBooleanValue());
//				}
//		});
//    	snpcalling.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				if(snpcalling.getBooleanValue()){
//				}else{
//					outvariantsonly.setBooleanValue(false);
//				}
//				}
//		});
//
//
//    	outvariantsonly.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				if(outvariantsonly.getBooleanValue()){
//					snpcalling.setBooleanValue(true);
//				}
//				}
//		});
//
//    	constrainedcalling.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				if(!constrainedcalling.getStringValue().equals("pair")&&!constrainedcalling.getStringValue().equals("No constrains")){
//					ifsamplelist.setBooleanValue(true);
//				}
//				}
//		}); 	
 
    	
    }
}

