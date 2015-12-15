package de.helmholtz_muenchen.ibis.ngs.geneSetAnalysis;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.CaseControlAnalyzerNodeDialog;

/**
 * <code>NodeDialog</code> for the "GeneSetAnalysis" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GeneSetAnalysisNodeDialog extends CaseControlAnalyzerNodeDialog {

    /**
     * New pane for configuring the GeneSetAnalysis node.
     */
    protected GeneSetAnalysisNodeDialog() {

    }

	@Override
	protected void addDialogComponent() {

	}
}

