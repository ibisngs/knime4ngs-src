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
		final SettingsModelString m_resolution = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_RESOLUTION, GeneBasedAnalysisNodeModel.RESOLUTION[0]);
		final SettingsModelBoolean m_fisher = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_FISHER,true);
	    final SettingsModelString m_alt_fisher = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ALT_FISHER,"two.sided");
	    final SettingsModelBoolean m_wilcoxon = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_WILCOXON,true);
	    final SettingsModelString m_alt_wilcoxon = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ALT_WILCOXON,"two.sided");
//		final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
//	    final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(GeneBasedAnalysisNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
	    final SettingsModelBoolean m_fisher_bg = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_FISHER_BACKGROUND,true);
	    final SettingsModelBoolean m_hyper = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_HYPER_BACKGROUND,true);
	    final SettingsModelString m_order_by = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ORDER_BY,GeneBasedAnalysisNodeModel.METHODS[0]);
		
	    createNewGroup("Statistic settings");
	    addDialogComponent(new DialogComponentStringSelection(m_resolution, "Choose resolution",GeneBasedAnalysisNodeModel.RESOLUTION));
	    addDialogComponent(new DialogComponentBoolean(m_fisher,"Perform Fisher's exact test?"));
	    addDialogComponent(new DialogComponentStringSelection(m_alt_fisher,"Choose alternative hypothesis", GeneBasedAnalysisNodeModel.ALTERNATIVES));
	    addDialogComponent(new DialogComponentBoolean(m_wilcoxon,"Compute Wilcoxon p-value?"));
	    addDialogComponent(new DialogComponentStringSelection(m_alt_wilcoxon,"Choose alternative hypothesis", GeneBasedAnalysisNodeModel.ALTERNATIVES));
//	    addDialogComponent(new DialogComponentBoolean(m_bin_back,"Compute binomial background?"));
//	    addDialogComponent(new DialogComponentNumber(m_pseudo_freq,"Background pseudo frequency",0.001));
	    addDialogComponent(new DialogComponentBoolean(m_fisher_bg, "Compute Fisher's exact versus background?"));
	    addDialogComponent(new DialogComponentBoolean(m_hyper,"Compute hypergeometric background?"));
	    addDialogComponent(new DialogComponentStringSelection(m_order_by,"Order by p-values of",GeneBasedAnalysisNodeModel.METHODS));
	    
	}
}

