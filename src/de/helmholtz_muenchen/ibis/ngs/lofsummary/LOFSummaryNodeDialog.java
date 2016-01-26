package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
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
	private final SettingsModelString cdsin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_CDS_INFILE,"-");
	private final SettingsModelString pedin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_PED_INFILE,"-");
	private final SettingsModelBoolean internal_gene_set = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_INTERNAL_GENE_SET,true);
	private final SettingsModelBoolean create_var_sum = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_CREATE_VAR_SUM,true);
//	private final SettingsModelBoolean create_gene_sum = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_CREATE_GENE_SUM,true);
	private final SettingsModelBoolean create_sample_sum = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_CREATE_SAMPLE_SUM,true);
	private final SettingsModelBoolean parallel_exec = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_PARALLEL_EXEC,false);
	private final SettingsModelBoolean create_matrix = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_CREATE_MATRIX,true);

	
    protected LOFSummaryNodeDialog() {
    	
    	createNewGroup("Path to CDS/GTF file");
    	addDialogComponent(new DialogComponentFileChooser(cdsin, "his_id_LOFStatistics_CDSIN", 0, ".fa|.fasta|.gtf"));
    	addDialogComponent(new DialogComponentBoolean(internal_gene_set, "Use internal representation of CDS/GTF file?"));
    	
    	createNewGroup("Path to PED file");
    	addDialogComponent(new DialogComponentFileChooser(pedin, "his_id_LOFStatistics_PEDIN", 0, ".ped"));
    	
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(create_var_sum,"Create variant summary?"));
//    	addDialogComponent(new DialogComponentBoolean(create_gene_sum,"Create gene summary?"));
    	addDialogComponent(new DialogComponentBoolean(create_sample_sum,"Create sample summary?"));
    	addDialogComponent(new DialogComponentBoolean(create_matrix, "Create matrix file?"));
    	addDialogComponent(new DialogComponentBoolean(parallel_exec,"Create summaries in parallel?"));
    }
}

