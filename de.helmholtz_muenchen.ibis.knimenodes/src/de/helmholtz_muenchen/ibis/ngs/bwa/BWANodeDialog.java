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

package de.helmholtz_muenchen.ibis.ngs.bwa;

import java.io.File;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;

/**
 * <code>NodeDialog</code> for the "BWA" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiterr
 * 
 */
public class BWANodeDialog extends HTExecutorNodeDialog {

	
	
    /**
     * New pane for configuring the BWA node.
     */
	protected BWANodeDialog() {
		super(FastQCell.TYPE.getPreferredValueClass(), FastQCell.TYPE.getPreferredValueClass());
	}
	
	public void addToolDialogComponents() {
		
		final SettingsModelString bwa 						= new SettingsModelString(BWANodeModel.CFGKEY_BWAFILE,"bwa");	
		final SettingsModelString refseq 					= new SettingsModelString(BWANodeModel.CFGKEY_REFSEQFILE,"");
		final SettingsModelBoolean indexrefseq 				= new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKINDEX, true);
		final SettingsModelString readGroup 				= new SettingsModelString(BWANodeModel.CFGKEY_READGROUP, "@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA");
		final SettingsModelBoolean readGroupBoolean 		= new SettingsModelBoolean(BWANodeModel.CFGKEY_READGROUPBOOLEAN, false);
		final SettingsModelIntegerBounded ALN_THREADS 		= new SettingsModelIntegerBounded(BWANodeModel.CFGKEY_THREADS,2, 1, Integer.MAX_VALUE);
//		final SettingsModelOptionalString Optional_Index 	= new SettingsModelOptionalString(BWANodeModel.CFGKEY_OPTIONAL_Index,"",false);
		final SettingsModelOptionalString Optional_Aln 		= new SettingsModelOptionalString(BWANodeModel.CFGKEY_OPTIONAL_Aln,"",false);
		final SettingsModelOptionalString Optional_Map 		= new SettingsModelOptionalString(BWANodeModel.CFGKEY_OPTIONAL_Map,"",false);
		final SettingsModelString alnalgo 					= new SettingsModelString(BWANodeModel.CFGKEY_ALNALGO,"BWA-MEM");
		final SettingsModelString bwtIndex 					= new SettingsModelString(BWANodeModel.CFGKEY_BWTINDEX,"BWT-SW");
//		final SettingsModelBoolean checkColorSpaced			= new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKCOLORSPACED, false);
		
    	addPrefPageSetting(bwa, IBISKNIMENodesPlugin.BWA);
    	addPrefPageSetting(refseq, IBISKNIMENodesPlugin.REF_GENOME);
    	readGroup.setEnabled(false);
    	
    	bwa.setEnabled(false);
    	addDialogComponent(new DialogComponentBoolean(indexrefseq, "Index reference sequence (Has to be done if index does not exist yet)."));
    	addDialogComponent(new DialogComponentStringSelection(bwtIndex,"Algorithm for constructing BWT index:","BWT-SW","IS"));
    	addDialogComponent(new DialogComponentStringSelection(alnalgo,"Algorithm for mapping:","BWA-MEM","BWA-backtrack","BWA-SW"));
//    	addDialogComponent(new DialogComponentBoolean(checkColorSpaced, "Build color-space index."));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(readGroupBoolean,"Specify read group header:"));
    	addDialogComponent(new DialogComponentString(readGroup,""));
    	setHorizontalPlacement(false);
//    	addDialogComponent(new DialogComponentOptionalString(Optional_Index, "Optional Indexing Parameters"));
    	addDialogComponent(new DialogComponentOptionalString(Optional_Aln, "Optional 'Aln' Parameters"));
    	Optional_Aln.setEnabled(false);
    	addDialogComponent(new DialogComponentOptionalString(Optional_Map, "Optional Mapping Parameters"));
    	addDialogComponent(new DialogComponentNumber(ALN_THREADS, "Number of threads", 1));
    	
    	refseq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
    	    	String path2refSeq = refseq.getStringValue();
    			File file0 = new File(path2refSeq+".amb");
    			File file1 = new File(path2refSeq+".ann");
    			File file2 = new File(path2refSeq+".bwt");
    			File file3 = new File(path2refSeq+".pac");
    			File file4 = new File(path2refSeq+".sa");
    			if(file0.exists() && file1.exists() && file2.exists() && file3.exists() && file4.exists()){
    				indexrefseq.setBooleanValue(false);
    				indexrefseq.setEnabled(true);
    			} else {
    				indexrefseq.setBooleanValue(true);
    				indexrefseq.setEnabled(false);
    			}
			}
		});
    	
    	readGroupBoolean.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(readGroupBoolean.getBooleanValue()){
						readGroup.setEnabled(true);
					}else{
						readGroup.setEnabled(false);
					}
			}
		});
    	
    	alnalgo.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(alnalgo.getStringValue().equals("BWA-backtrack")){
						Optional_Aln.setEnabled(true);
					}else{
						Optional_Aln.setEnabled(false);
					}
			}
		});
    	
    	
    	
    }	
}

