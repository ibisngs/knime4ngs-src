package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "CaseControlAnalyzer" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class CaseControlAnalyzerNodeDialog extends DefaultNodeSettingsPane {

	//model file settings models
    private final SettingsModelString m_model_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_MODEL_GENE_ID, "gene_id");
    private final SettingsModelString m_freq = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_FREQ, "lof_freq");
    private final SettingsModelDouble m_pseudo_freq = new SettingsModelDouble(CaseControlAnalyzerNodeModel.CFGKEY_PSEUDO_FREQ,0.0);
    
    //summary file settings models
    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_ctrl");
    /**
     * New pane for configuring the CaseControlAnalyzer node.
     */
    protected CaseControlAnalyzerNodeDialog() {
    	
    	createNewGroup("Summary file settings");
    	addDialogComponent(new DialogComponentString(m_summary_gene_id,"Gene identifier column"));
    	addDialogComponent(new DialogComponentString(m_case_cond, "Case condition column"));
    	addDialogComponent(new DialogComponentString(m_case_ncond, "Case non-condition column"));
    	addDialogComponent(new DialogComponentString(m_control_cond, "Control condition column"));
    	addDialogComponent(new DialogComponentString(m_control_ncond, "Control non-condition column"));
    	
    	createNewGroup("Model (file) settings");
    	addDialogComponent(new DialogComponentString(m_model_gene_id, "Model gene identifier column"));
    	addDialogComponent(new DialogComponentString(m_freq, "Gene-based frequency column"));
    	addDialogComponent(new DialogComponentNumber(m_pseudo_freq,"Pseudo frequency",0.001));
    	
    	//TODO add fields to choose which statistic(s) to be computed and output ordered by which statistic

    }
}

