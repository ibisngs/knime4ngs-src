package de.helmholtz_muenchen.ibis.ngs.geneSetAnalysis;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

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

    	final SettingsModelString genesetin = new SettingsModelString(GeneSetAnalysisNodeModel.CFGKEY_GENE_SET_INFILE,"");
    	final SettingsModelString genesetname = new SettingsModelString(GeneSetAnalysisNodeModel.CFGKEY_GENE_SET,"gene_set");

    	
		createNewGroup("Gene set file");
    	addDialogComponent(new DialogComponentFileChooser(genesetin,"his_id_LOFStatistics_GENESET",0,".gmt"));
    	addDialogComponent(new DialogComponentString(genesetname,"Gene set symbol"));
	}
}

