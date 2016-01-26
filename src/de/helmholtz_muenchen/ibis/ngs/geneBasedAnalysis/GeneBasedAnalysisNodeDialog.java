package de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.CaseControlAnalyzerNodeDialog;

/**
 * <code>NodeDialog</code> for the "GeneBasedAnalysis" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GeneBasedAnalysisNodeDialog extends CaseControlAnalyzerNodeDialog {

    /**
     * New pane for configuring the GeneBasedAnalysis node.
     */
    protected GeneBasedAnalysisNodeDialog() {

    }

	@Override
	protected void addDialogComponent() {
		final SettingsModelBoolean m_fisher = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_FISHER,true);
	    final SettingsModelString m_alternative = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ALT,GeneBasedAnalysisNodeModel.ALTERNATIVES[0]);
//	    final SettingsModelBoolean m_wilcoxon = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_WILCOXON,true);
//		final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
//	    final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(GeneBasedAnalysisNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
	    final SettingsModelBoolean m_hyper = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_HYPER_BACKGROUND,true);
	    final SettingsModelString m_order_by = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ORDER_BY,GeneBasedAnalysisNodeModel.METHODS[0]);
		
	    createNewGroup("Statistic settings");
	    addDialogComponent(new DialogComponentBoolean(m_fisher,"Perform Fisher's exact test?"));
	    addDialogComponent(new DialogComponentStringSelection(m_alternative,"Choose alternative hypothesis", GeneBasedAnalysisNodeModel.ALTERNATIVES));
//	    addDialogComponent(new DialogComponentBoolean(m_wilcoxon,"Compute Wilcoxon p-value?"));
//	    addDialogComponent(new DialogComponentBoolean(m_bin_back,"Compute binomial background?"));
//	    addDialogComponent(new DialogComponentNumber(m_pseudo_freq,"Background pseudo frequency",0.001));
	    addDialogComponent(new DialogComponentBoolean(m_hyper,"Compute hypergeometric background?"));
	    addDialogComponent(new DialogComponentStringSelection(m_order_by,"Order by p-values of",GeneBasedAnalysisNodeModel.METHODS));
	    
	}
}

