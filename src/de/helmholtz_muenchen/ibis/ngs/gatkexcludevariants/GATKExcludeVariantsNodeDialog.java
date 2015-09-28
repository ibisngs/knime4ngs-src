package de.helmholtz_muenchen.ibis.ngs.gatkexcludevariants;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKExcludeVariants" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GATKExcludeVariantsNodeDialog extends SelectVariantsNodeDialog {


	@Override
	protected void addFilterDialogComponent() {
		addDialogComponent(new DialogComponentString(getFILTERSTRINGModel(), "Parameters"));
	}
}

