package de.helmholtz_muenchen.ibis.ngs.gatkgenotypeconcordance;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKGenotypeConcordance" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilan Hastreiter
 */
public class GATKGenotypeConcordanceNodeDialog extends GATKNodeDialog {
	
    @Override
	protected void addDialogComponent() {
    	
        final SettingsModelString EVAL = new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_EVAL, "");
        final SettingsModelString COMP = new SettingsModelString(GATKGenotypeConcordanceNodeModel.CFGKEY_COMP, "");
    	final SettingsModelBoolean UNSAFE = new SettingsModelBoolean(GATKGenotypeConcordanceNodeModel.CFGKEY_UNSAFE,false);

        
    	createNewGroup("The variants and genotypes to evaluate");
    	addDialogComponent(new DialogComponentFileChooser(EVAL, "eval", JFileChooser.OPEN_DIALOG, false, ".vcf"));
 	
    	createNewGroup("The variants and genotypes to compare against");
    	addDialogComponent(new DialogComponentFileChooser(COMP, "comp", JFileChooser.OPEN_DIALOG, false, ".vcf"));
    	
    	addDialogComponent(new DialogComponentBoolean(UNSAFE, "Enable unsafe-mode (not recommended!)"));
		
	}
}

