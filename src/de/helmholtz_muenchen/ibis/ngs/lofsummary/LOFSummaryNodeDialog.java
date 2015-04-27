package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "LOFStatistics" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class LOFSummaryNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the LOFStatistics node.
     */
	private final SettingsModelString vcfin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_VCF_INFILE,"-");
	private final SettingsModelString cdsin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_CDS_INFILE,"-");
	private final SettingsModelString annotation = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_ANNOTATION, "");
	private final SettingsModelBoolean high_confidence = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_HIGH_CONFIDENCE,false);


    protected LOFSummaryNodeDialog() {
    	
    	//input file
    	createNewGroup("Path to VCF file");
    	addDialogComponent(new DialogComponentFileChooser(vcfin, "his_id_LOFStatistics_VCFIN", 0, ".vcf"));
    	
    	createNewGroup("Path to CDS file");
    	addDialogComponent(new DialogComponentFileChooser(cdsin, "his_id_LOFStatistics_CDSIN", 0, ".fa"));
    	
    	//annotation selection
        createNewGroup("Annotation Selection");
        addDialogComponent(new DialogComponentStringSelection(annotation, "Tool", LOFSummaryNodeModel.ANNOTATIONS_AVAILABLE));
        addDialogComponent(new DialogComponentBoolean(high_confidence, "Restrict to high confidence LoFs?"));
    }
}

