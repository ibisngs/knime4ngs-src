package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GeneticBackgroundModel" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeDialog extends DefaultNodeSettingsPane {
	
	private final SettingsModelString gtf_aff = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_GTF_AFF,GeneticBackgroundModelNodeModel.BASIS[0]);
	private final SettingsModelString m_ac = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AC,"AC");
	private final SettingsModelString m_an = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AN,"AN");
	private final SettingsModelString genesetin = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_GENE_SET_INFILE,"");
    private final SettingsModelBoolean use_symbol = new SettingsModelBoolean(GeneticBackgroundModelNodeModel.CFGKEY_USE_SYMBOL, false);

	
    /**
     * New pane for configuring the GeneticBackgroundModel node.
     */
    protected GeneticBackgroundModelNodeDialog() {
    	
    	createNewGroup("Compute gene-(set)-based variant frequencies using");
    	addDialogComponent(new DialogComponentButtonGroup(gtf_aff, true, null, GeneticBackgroundModelNodeModel.BASIS));
    	
    	createNewGroup("Field IDs indicating the AC and AN to be used");
    	m_ac.setEnabled(false);
    	m_an.setEnabled(false);
    	addDialogComponent(new DialogComponentString(m_ac, "Allele Count (AC)"));
    	addDialogComponent(new DialogComponentString(m_an, "Allele Number (AN)"));
    	
    	gtf_aff.addChangeListener(new ChangeListener () {

			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(gtf_aff.getStringValue().equals(GeneticBackgroundModelNodeModel.BASIS[0])) {
					m_ac.setEnabled(false);
			    	m_an.setEnabled(false);
				} else {
					m_ac.setEnabled(true);
			    	m_an.setEnabled(true);
				}
			}
    	});
    	
    	createNewGroup("Path to gene set file");
    	addDialogComponent(new DialogComponentFileChooser(genesetin,"his_id_gb_GENESET",0,".gmt"));
    	
    	createNewGroup("Output options");
    	addDialogComponent(new DialogComponentBoolean(use_symbol, "Use gene symbols as unique identifier?"));
    }
}

