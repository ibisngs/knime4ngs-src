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


package de.helmholtz_muenchen.ibis.ngs.fastSam2Bam;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * Dialog for a faster implementation of a sam to bam converter.
 * @author Michael Kluge
 */
public class FastSam2BamNodeDialog extends HTExecutorNodeDialog {

	
	
    
    /**
     * New pane for configuring the FastSam2Bam node.
     */
    protected FastSam2BamNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	// definition of SettingsModel (all prefixed with SET)
    	final SettingsModelString SET_GENOME 			= new SettingsModelString(FastSam2BamNodeModel.CFGKEY_GENOME_FILE, FastSam2BamNodeModel.DEFAULT_GENOME_FILE);
        final SettingsModelString SET_OUTPUT_PATH		= new SettingsModelString(FastSam2BamNodeModel.CFGKEY_OUTPUT_PATH, FastSam2BamNodeModel.DEFAULT_OUTPUT_PATH);
        final SettingsModelInteger SET_CORE_NUMBER		= new SettingsModelInteger(FastSam2BamNodeModel.CFGKEY_CORE_NUMBER, FastSam2BamNodeModel.DEFAULT_CORE_NUMBER);
        final SettingsModelInteger SET_SPLIT_SIZE		= new SettingsModelInteger(FastSam2BamNodeModel.CFGKEY_SPLIT_SIZE, FastSam2BamNodeModel.DEFAULT_SPLIT_SIZE);
        final SettingsModelString SET_TMP_PATH			= new SettingsModelString(FastSam2BamNodeModel.CFGKEY_TMP_PATH, FastSam2BamNodeModel.DEFAULT_TMP_PATH);
        final SettingsModelBoolean SET_USE_RAM_AS_TMP  	= new SettingsModelBoolean(FastSam2BamNodeModel.CFGKEY_USE_RAM, FastSam2BamNodeModel.DEFAULT_USE_RAM);
        final SettingsModelString SET_PATH_SAMTOOLS		= new SettingsModelString(FastSam2BamNodeModel.CFGKEY_PATH_SAMTOOLS, FastSam2BamNodeModel.DEFAULT_PATH_SAMTOOLS);
//        private final SettingsModelString SET_PATH_PICTOOLS    	= new SettingsModelString(FastSam2BamNodeModel.CFGKEY_PATH_PICTOOLS, FastSam2BamNodeModel.DEFAULT_PATH_PICTOOLS);
        final SettingsModelBoolean SET_DELETE_SAM		= new SettingsModelBoolean(FastSam2BamNodeModel.CFGKEY_DELETE_SAM, FastSam2BamNodeModel.DEFAULT_DELETE_SAM);
    	
        // create open file/folder components
//        DialogComponentFileChooser dcPathSamtools 	= new DialogComponentFileChooser(SET_PATH_SAMTOOLS, "his_id_fs2b_samtools", 0, false);
//        DialogComponentFileChooser dcPathPictools 	= new DialogComponentFileChooser(SET_PATH_PICTOOLS, "his_id_fs2b_pictools", 0, true);
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_PATH, "his_id_fs2b_outpath", 0, true);
//        DialogComponentFileChooser dcGenomeFile 	= new DialogComponentFileChooser(SET_GENOME, "his_id_fs2b_genomeFile", 0, false);
       	DialogComponentString dcTmpPath				= new DialogComponentString(SET_TMP_PATH, "temp path:"); 
       	DialogComponentNumber dcThreadNumber		= new DialogComponentNumber(SET_CORE_NUMBER, "number of cores", 1);
       	DialogComponentNumber dcSplitSize			= new DialogComponentNumber(SET_SPLIT_SIZE, "number of reads per split", 1);
     	DialogComponentBoolean dcUseRamAsTmp 		= new DialogComponentBoolean(SET_USE_RAM_AS_TMP, "use RAM (" + FastSam2BamNodeModel.DEFAULT_USE_RAM_PATH + ") as temp folder");
     	DialogComponentBoolean dcDeleteSam			= new DialogComponentBoolean(SET_DELETE_SAM, "delete sam file after conversion");
      
    	addPrefPageSetting(SET_PATH_SAMTOOLS, IBISKNIMENodesPlugin.SAMTOOLS);
    	addPrefPageSetting(SET_GENOME, IBISKNIMENodesPlugin.REF_GENOME);
     	
       	// set a new title to them
     	dcOutputFolder.setBorderTitle("path to output folder");
//     	dcPathSamtools.setBorderTitle("path to samtool binaries");
//     	dcPathPictools.setBorderTitle("path to picard tool binaries");
//     	dcGenomeFile.setBorderTitle("path to genome in fasta format");
     
       	// add groups and components
//        createNewGroup("path to binaries");
//        addDialogComponent(dcPathSamtools);
//        addDialogComponent(dcPathPictools);
        
        createNewGroup("output");
        addDialogComponent(dcOutputFolder);
        
        createNewGroup("further options");
//        addDialogComponent(dcGenomeFile);
        addDialogComponent(dcThreadNumber);
        addDialogComponent(dcSplitSize);
        addDialogComponent(dcTmpPath);
        addDialogComponent(dcUseRamAsTmp);
        
        createNewGroup("WARNING: DO USE WITH CARE!");
        addDialogComponent(dcDeleteSam);
        
        // some event listeners
        SET_USE_RAM_AS_TMP.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(SET_USE_RAM_AS_TMP.getBooleanValue()) 
					SET_TMP_PATH.setStringValue(FastSam2BamNodeModel.DEFAULT_USE_RAM_PATH);
			}	
        });
        
        SET_CORE_NUMBER.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(SET_CORE_NUMBER.getIntValue() < 2)
					SET_CORE_NUMBER.setIntValue(2);
			}	
        });
        
        SET_SPLIT_SIZE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				//if(SET_SPLIT_SIZE.getIntValue() < 50000)
				//	SET_SPLIT_SIZE.setIntValue(50000);
			}	
        });
    }
    
//	protected void updatePrefs() {
//		if(usePrefPage.getBooleanValue()) {
//	    	String toolPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("samtools");
//	    	if(toolPath != null && !toolPath.equals("")) {
//	    		SET_PATH_SAMTOOLS.setStringValue(toolPath);
//	    		SET_PATH_SAMTOOLS.setEnabled(false);
//	    	} else {
//	    		SET_PATH_SAMTOOLS.setEnabled(true);
//	    	}
//	    	
//		} else {
//			SET_PATH_SAMTOOLS.setEnabled(true);
//		}
//	}
    
//    public void onOpen() {
//    	String toolPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("samtools");
//    	if(toolPath != null && !toolPath.equals("")) {
//    		SET_PATH_SAMTOOLS.setStringValue(toolPath);
//    	}
//    	String refGenome = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
//    	if(refGenome != null && !refGenome.equals("")) {
//    		SET_GENOME.setStringValue(refGenome);
//    	}
//    }
}

