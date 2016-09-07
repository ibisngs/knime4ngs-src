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
package de.helmholtz_muenchen.ibis.ngs.featureCounts;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "FeatureCounts" Node.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeDialog extends HTExecutorNodeDialog {

	private final static String BINARY_NAME = IBISKNIMENodesPlugin.FEATURE_COUNTS;
	
	

    
    protected FeatureCountsNodeDialog() {}
    
    public void addToolDialogComponents() {
        
    	// definition of SettingsModel (all prefixed with SET)
        final SettingsModelString SET_FEATURE_TYPE				= new SettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FEATURE, FeatureCountsNodeModel.DEFAULT_ANNOTATION_FEATURE);
        final SettingsModelString SET_OUTPUT_FILE 				= new SettingsModelString(FeatureCountsNodeModel.CFGKEY_OUTPUT_FOLDER, "");
        final SettingsModelString SET_ANNOTATION_FILE			= new SettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FILE, FeatureCountsNodeModel.DEFAULT_ANNOTATION_FILE);
        final SettingsModelString SET_ANNOTATION_TYPE			= new SettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_TYPE, FeatureCountsNodeModel.DEFAULT_ANNOTATION_TYPE);
        final SettingsModelInteger SET_THREAD_NUMBER			= new SettingsModelInteger(FeatureCountsNodeModel.CFGKEY_THREAD_NUMBER, FeatureCountsNodeModel.DEFAULT_THREAD_NUMBER);
        final SettingsModelBoolean SET_COUNT_MULTIMAPPED		= new SettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_MULTIMAPPED, FeatureCountsNodeModel.DEAFULT_COUNT_MULTIMAPPED);
        final SettingsModelBoolean SET_COUNT_OVERLAPPING_MULTI	= new SettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_OVERLAPPING_MULTI, FeatureCountsNodeModel.DEAFULT_COUNT_MULTI_OVERLAPING);
        final SettingsModelBoolean SET_COUNT_FRAGMENTS			= new SettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_FRAGMENTS, FeatureCountsNodeModel.DEAFULT_COUNT_FRAGMENTS);
        final SettingsModelBoolean SET_CHIMERIC_FRAGMENTS		= new SettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_CHIMERIC_FRAGMENTS, FeatureCountsNodeModel.DEAFULT_COUNT_CHIMERIC_FRAGMENTS);
        final SettingsModelBoolean SET_FEATURE_LEVEL			= new SettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_ON_FEATURE_LVL, FeatureCountsNodeModel.DEFAULT_COUNT_ON_FEATURE_LVL);
        final SettingsModelString SET_GROUP_FEATURE				= new SettingsModelString(FeatureCountsNodeModel.CFGKEY_GROUP_FEATURE, FeatureCountsNodeModel.DEFAULT_GROUP_FEATURE);
        final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);

        
        // create open file/folder components
        DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE_FeatureCounts", 0, true);
       	DialogComponentFileChooser dcAnnotationFile = new DialogComponentFileChooser(SET_ANNOTATION_FILE, "his_id_GENOME_FOLDER_FeatureCounts", 0, ".gtf", ".saf");
       	DialogComponentNumber dcThreadNumber		= new DialogComponentNumber(SET_THREAD_NUMBER, "thread number", 1);
       	DialogComponentString dcFeatureType			= new DialogComponentString(SET_FEATURE_TYPE, "feature type used for counting:"); 					
     	DialogComponentBoolean dcCountMultimapped 	= new DialogComponentBoolean(SET_COUNT_MULTIMAPPED, "count multimapped reads");
     	DialogComponentBoolean dcCountOverlapping 	= new DialogComponentBoolean(SET_COUNT_OVERLAPPING_MULTI, "count reads that overlapp multiple features");
     	DialogComponentBoolean dcCountFragments 	= new DialogComponentBoolean(SET_COUNT_FRAGMENTS, "count fragments instead of reads");
     	DialogComponentBoolean dcCountChimeric	 	= new DialogComponentBoolean(SET_CHIMERIC_FRAGMENTS, "count chimeric (paired reads)");
     	DialogComponentBoolean dcCountOnFeatureLvl	= new DialogComponentBoolean(SET_FEATURE_LEVEL, "perform read summarization at the feature level (eg. exon level) ");
       	DialogComponentString dcGroupFeature		= new DialogComponentString(SET_GROUP_FEATURE, "feature type used for grouping results:"); 
     	
       	// create string selection component
       	DialogComponentStringSelection dcAnnotationType 	= new DialogComponentStringSelection(SET_ANNOTATION_TYPE, "file type:", FeatureCountsNodeModel.DEFAULT_ANNOTATION_TYPE, FeatureCountsNodeModel.ALTERNATIVE_ANNOTATION_TYPE);
       	
        addPrefPageSetting(SET_BINARY_PATH,getNameOfBinary());
//       	createNewTab("FeatureCounts");
       	// set a new title to them
       	dcOutputFile.setBorderTitle("Path to output folder");
       	dcAnnotationFile.setBorderTitle("Path to annotation file");
     
       	// add groups and components
        createNewGroup("Input");
        addDialogComponent(dcAnnotationFile);
        
        createNewGroup("Output");
        addDialogComponent(dcOutputFile);
        
        createNewGroup("Further options");
        addDialogComponent(dcAnnotationType);
        addDialogComponent(dcFeatureType);
        addDialogComponent(dcCountMultimapped);
        addDialogComponent(dcCountOverlapping);
        addDialogComponent(dcThreadNumber);
        addDialogComponent(dcCountOnFeatureLvl);
        addDialogComponent(dcGroupFeature);
        
        createNewGroup("Paired read options");
        addDialogComponent(dcCountFragments);
        addDialogComponent(dcCountChimeric);
        
        // 
        SET_THREAD_NUMBER.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				// check range
				int threads = SET_THREAD_NUMBER.getIntValue();
				if(FeatureCountsNodeModel.MIN_THREADS > threads)
					SET_THREAD_NUMBER.setIntValue(FeatureCountsNodeModel.MIN_THREADS);
				if(FeatureCountsNodeModel.MAX_THREADS < threads)
					SET_THREAD_NUMBER.setIntValue(FeatureCountsNodeModel.MAX_THREADS);
			}	
        });
    }


	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}

