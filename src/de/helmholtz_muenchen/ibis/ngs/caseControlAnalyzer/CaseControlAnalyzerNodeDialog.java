package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
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
    private final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(CaseControlAnalyzerNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
    
    //summary file settings models
    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_NCOND, "un_ctrl");
    
    //method settings models
    private final SettingsModelBoolean m_fisher = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_FISHER,true);
    private final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
    private final SettingsModelString m_order_by = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_ORDER_BY,CaseControlAnalyzerNodeModel.METHODS[0]);
    
    
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
    	
    	createNewGroup("Choose statistics and ordering");
    	addDialogComponent(new DialogComponentBoolean(m_fisher, "Compute Fisher's Exact?"));
    	addDialogComponent(new DialogComponentBoolean(m_bin_back, "Compute Binomial Background?"));
    	addDialogComponent(new DialogComponentStringSelection(m_order_by, "Order by", CaseControlAnalyzerNodeModel.METHODS));
    }
}

