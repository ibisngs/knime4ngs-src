package de.helmholtz_muenchen.ibis.ngs.fastqc_v2;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.Bcl2FastQ.Bcl2FastQNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "FastQC_v2" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Paul Hager
 */
public class FastQC_v2NodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the FastQC_v2 node.
     */
    protected FastQC_v2NodeDialog() {

    }

	@Override
	public void addToolDialogComponents() {
		//options tab
		final SettingsModelString outfolder = new SettingsModelString(FastQC_v2NodeModel.CFGKEY_OUTFOLDER,"");
		final SettingsModelString fastqc = new SettingsModelString(FastQC_v2NodeModel.CFGKEY_FASTQC,"");
		final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(FastQC_v2NodeModel.CFGKEY_THREADS,4, 1, Integer.MAX_VALUE);
		
		createNewGroup("Folder for output files");
		addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_FastQC_v2_OUT", 0, true));
		
		createNewGroup("Number of threads");
		addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1));
		
		addPrefPageSetting(fastqc, IBISKNIMENodesPlugin.FASTQC);
	}
}

