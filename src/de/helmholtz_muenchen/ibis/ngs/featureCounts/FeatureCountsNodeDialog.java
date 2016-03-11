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
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeDialog;

/**
 * <code>NodeDialog</code> for the "FeatureCounts" Node.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeDialog extends BinaryWrapperNodeDialog {

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

	@Override
	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}

