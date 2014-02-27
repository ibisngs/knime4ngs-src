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
 * <code>NodeDialog</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeDialog extends BinaryWrapperNodeDialog {

	private final static String BINARY_NAME = "featureCounts";
	
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FEATURE_TYPE			= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FEATURE, this);
    private final SettingsModelString SET_OUTPUT_FILE 			= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_OUTPUT_FILE, this);
    private final SettingsModelString SET_ANNOTATION_FILE		= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_FILE, this);
    private final SettingsModelString SET_ANNOTATION_TYPE		= FeatureCountsNodeModel.getSettingsModelString(FeatureCountsNodeModel.CFGKEY_ANNOTATION_TYPE, this);
    private final SettingsModelInteger SET_THREAD_NUMBER		= FeatureCountsNodeModel.getSettingsModelInteger(FeatureCountsNodeModel.CFGKEY_THREAD_NUMBER, this);
    private final SettingsModelBoolean SET_COUNT_MULTIMAPPED	= FeatureCountsNodeModel.getSettingsModelBoolean(FeatureCountsNodeModel.CFGKEY_COUNT_MULTIMAPPED, this);
    
    protected FeatureCountsNodeDialog() {
        super();
       
        // create open file/folder components
        DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE", 0, true);
       	DialogComponentFileChooser dcAnnotationFile = new DialogComponentFileChooser(SET_ANNOTATION_FILE, "his_id_GENOME_FOLDER", 0, true);
       	DialogComponentNumber dcThreadNumber		= new DialogComponentNumber(SET_THREAD_NUMBER, "thread number", 1);
       	DialogComponentString dcFeatureType			= new DialogComponentString(SET_FEATURE_TYPE, "feature type used for counting:"); 					
       	DialogComponentBoolean dcCountMultimapped 	= new DialogComponentBoolean(SET_COUNT_MULTIMAPPED, "count multimapped reads");
       			
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
        addDialogComponent(dcThreadNumber);
        
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

