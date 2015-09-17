package de.helmholtz_muenchen.ibis.ngs.gatkfiltervqslod;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKFilterVQSLOD" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKFilterVQSLODNodeDialog extends SelectVariantsNodeDialog {

	@Override
	protected void addFilterDialogComponent() {
		addDialogComponent(new DialogComponentString(getFILTERSTRINGModel(), "Enter VQSLOD Cutoff"));

		
	}
	
}

