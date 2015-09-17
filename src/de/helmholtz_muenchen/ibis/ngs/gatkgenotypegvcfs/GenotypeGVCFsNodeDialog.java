package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GenotypeGVCFs" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeDialog extends GATKNodeDialog {

	
	

	
	@Override
	protected void addDialogComponent() {
		
		final SettingsModelIntegerBounded NT = new SettingsModelIntegerBounded(GenotypeGVCFsNodeModel.CFGKEY_NT_FILE, 1, 1, Integer.MAX_VALUE);
		
		
		addDialogComponent(new DialogComponentNumber(NT, "Number of threads", 1));
		
		
	}

}

