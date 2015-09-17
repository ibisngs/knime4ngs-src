package de.helmholtz_muenchen.ibis.ngs.gatkselectvariants;


import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeDialog;
/**
 * <code>NodeDialog</code> for the "GATKSelectVariants" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class GATKSelectVariantsNodeDialog extends SelectVariantsNodeDialog {


	@Override
	protected void addFilterDialogComponent() {
		addDialogComponent(new DialogComponentStringSelection(getFILTERSTRINGModel(), "Select Variant Type", "SNP","INDEL"));
		
	}
}

