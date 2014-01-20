package de.helmholtz_muenchen.ibis.ngs.genscan;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Genscan" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class GenscanNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the Genscan node.
     */
    protected GenscanNodeDialog() {

    	final SettingsModelBoolean usematrix = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_USEMATRIX, true);
    	final SettingsModelString choosematrix = new SettingsModelString(GenscanNodeModel.CFGKEY_MATRIXFILE,null);
    	final SettingsModelBoolean subopt = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_SUBOPT, false);
    	final SettingsModelDoubleBounded suboptcutoff = new SettingsModelDoubleBounded(GenscanNodeModel.CFGKEY_SUBOPTCUTOFF, 0.1, 0.01, 0.99);
    	final SettingsModelBoolean ps = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_PS, true);
        final SettingsModelIntegerBounded psscale = new SettingsModelIntegerBounded(GenscanNodeModel.CFGKEY_PSSCALE, 1, 1, 100);
    	choosematrix.setEnabled(false);
    	suboptcutoff.setEnabled(false);
    	
    	createNewGroup("Genscan");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(GenscanNodeModel.CFGKEY_GENSCANFILE,""), "his_gene_id", 0, ""));
    	//createNewGroup("Sequence/ genome file (FastA)");
    	//addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(GenscanNodeModel.CFGKEY_SEQFILE,null), "his0_id", 0, null));
    	createNewGroup("Parameter matrix");
    	addDialogComponent(new DialogComponentBoolean(usematrix, "Use HumanIso matrix, located in the installpath of Genscan"));
    	addDialogComponent(new DialogComponentFileChooser(choosematrix, "his1_gene_id", 0, "smat"));
    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(GenscanNodeModel.CFGKEY_VERBOSEOUTPUT, true), "Add some extra explanatory information to the text output."));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(GenscanNodeModel.CFGKEY_CDS, true), "Also print predicted CDS."));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(subopt, "Identify suboptimal exons."));
    	addDialogComponent(new DialogComponentNumber(suboptcutoff, "Probability cutoff:", 0.01));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ps, "Create the PostScript (graphical) output."));
    	addDialogComponent(new DialogComponentNumber(psscale, "Scale to make the PostScript image:", 1));
    	setHorizontalPlacement(false);
    	
    	
    	usematrix.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				choosematrix.setEnabled(!usematrix.getBooleanValue());
			}
		});
    	
    	subopt.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				suboptcutoff.setEnabled(subopt.getBooleanValue());
			}
		});
    	
    	ps.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				psscale.setEnabled(ps.getBooleanValue());
			}
		});
    	
    }
}

