package de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
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
public abstract class CaseControlAnalyzerNodeDialog extends DefaultNodeSettingsPane {

	//model file settings models
    private final SettingsModelString m_model_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_MODEL_GENE_ID, "gene_id");
    private final SettingsModelString m_freq = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_FREQ, "variant_freq");
    private final SettingsModelIntegerBounded pop_size = new SettingsModelIntegerBounded(CaseControlAnalyzerNodeModel.CFGKEY_POP_SIZE, 1, 1, Integer.MAX_VALUE);
    
    //summary file settings models
//    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
//    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
//    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
//    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
//    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_NCOND, "un_ctrl");
//    private final SettingsModelBoolean m_ignore_unobserved = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_IGNORE_UNOBSERVED,false);
    
    /**
     * New pane for configuring the CaseControlAnalyzer node.
     */
    protected CaseControlAnalyzerNodeDialog() {
    	
//    	createNewGroup("Summary file settings");
//    	addDialogComponent(new DialogComponentString(m_summary_gene_id,"Gene identifier column"));
//    	addDialogComponent(new DialogComponentString(m_case_cond, "Case condition column"));
//    	addDialogComponent(new DialogComponentString(m_case_ncond, "Case non-condition column"));
//    	addDialogComponent(new DialogComponentString(m_control_cond, "Control condition column"));
//    	addDialogComponent(new DialogComponentString(m_control_ncond, "Control non-condition column"));
//    	addDialogComponent(new DialogComponentBoolean(m_ignore_unobserved, "Ignore genes/transcripts unaffected in cases and controls?"));
    	
    	createNewGroup("Statistic settings");
    	addDialogComponent();
    	
    	createNewGroup("Model file settings");
    	addDialogComponent(new DialogComponentString(m_model_gene_id, "Gene/transcript identifier column"));
    	addDialogComponent(new DialogComponentString(m_freq, "Background frequency column"));
	    addDialogComponent(new DialogComponentNumber(pop_size, "Background population size",100));
    }
    
    protected abstract void addDialogComponent();
}

