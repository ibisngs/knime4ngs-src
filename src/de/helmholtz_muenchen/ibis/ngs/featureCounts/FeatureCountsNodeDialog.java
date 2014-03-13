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

	private final static String BINARY_NAME = "featureCounts";
	
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FEATURE_TYPE				= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FEATURE, this);
    private final SettingsModelString SET_OUTPUT_FILE 				= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_OUTPUT_FILE, this);
    private final SettingsModelString SET_ANNOTATION_FILE			= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FILE, this);
    private final SettingsModelString SET_ANNOTATION_TYPE			= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_TYPE, this);
    private final SettingsModelInteger SET_THREAD_NUMBER			= FeatureCountsNodeModel.getSettingsModelInteger(FeatureCountsNodeModel.CFGKEY_THREAD_NUMBER, this);
    private final SettingsModelBoolean SET_COUNT_MULTIMAPPED		= FeatureCountsNodeModel.getSettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_MULTIMAPPED, this);
    private final SettingsModelBoolean SET_COUNT_OVERLAPPING_MULTI	= FeatureCountsNodeModel.getSettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_OVERLAPPING_MULTI, this);
    private final SettingsModelBoolean SET_COUNT_FRAGMENTS			= FeatureCountsNodeModel.getSettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_FRAGMENTS, this);
    private final SettingsModelBoolean SET_CHIMERIC_FRAGMENTS		= FeatureCountsNodeModel.getSettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_CHIMERIC_FRAGMENTS, this);
    
    protected FeatureCountsNodeDialog() {
        super();
       
        // create open file/folder components
        DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE_FeatureCounts", 0, false);
       	DialogComponentFileChooser dcAnnotationFile = new DialogComponentFileChooser(SET_ANNOTATION_FILE, "his_id_GENOME_FOLDER_FeatureCounts", 0, ".gtf", ".saf");
       	DialogComponentNumber dcThreadNumber		= new DialogComponentNumber(SET_THREAD_NUMBER, "thread number", 1);
       	DialogComponentString dcFeatureType			= new DialogComponentString(SET_FEATURE_TYPE, "feature type used for counting:"); 					
     	DialogComponentBoolean dcCountMultimapped 	= new DialogComponentBoolean(SET_COUNT_MULTIMAPPED, "count multimapped reads");
     	DialogComponentBoolean dcCountOverlapping 	= new DialogComponentBoolean(SET_COUNT_OVERLAPPING_MULTI, "count reads that overlapp multiple features");
     	DialogComponentBoolean dcCountFragments 	= new DialogComponentBoolean(SET_COUNT_FRAGMENTS, "count fragments instead of reads");
     	DialogComponentBoolean dcCountChimeric	 	= new DialogComponentBoolean(SET_CHIMERIC_FRAGMENTS, "count chimeric (paired reads)");
       	
       	// create string selection component
       	DialogComponentStringSelection dcAnnotationType 	= new DialogComponentStringSelection(SET_ANNOTATION_TYPE, "file type:", FeatureCountsNodeModel.DEFAULT_ANNOTATION_TYPE, FeatureCountsNodeModel.ALTERNATIVE_ANNOTATION_TYPE);
       	
       	// set a new title to them
       	dcOutputFile.setBorderTitle("path to output file");
       	dcAnnotationFile.setBorderTitle("path to annotation file");
     
       	// add groups and components
        createNewGroup("input");
        addDialogComponent(dcAnnotationFile);
        
        createNewGroup("output");
        addDialogComponent(dcOutputFile);
        
        createNewGroup("further options");
        addDialogComponent(dcAnnotationType);
        addDialogComponent(dcFeatureType);
        addDialogComponent(dcCountMultimapped);
        addDialogComponent(dcCountOverlapping);
        addDialogComponent(dcThreadNumber);
        
        createNewGroup("paired read options");
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

