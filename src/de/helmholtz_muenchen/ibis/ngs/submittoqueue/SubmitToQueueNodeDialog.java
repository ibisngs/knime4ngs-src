package de.helmholtz_muenchen.ibis.ngs.submittoqueue;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "SubmitToQueue" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class SubmitToQueueNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the SubmitToQueue node.
     */
    protected SubmitToQueueNodeDialog() {

    	
    	//TODO: qname optional angeben
    	addDialogComponent(new DialogComponentNumber(new SettingsModelIntegerBounded(SubmitToQueueNodeModel.CFGKEY_MEMORY, SubmitToQueueNodeModel.DEFAULT_MEMORY, 1, 100), 
    			"Reservate memory [GB]: ", /*step*/ 1));
    	addDialogComponent(new DialogComponentString(new SettingsModelString(SubmitToQueueNodeModel.CFGKEY_JOBNAME, "NGS"), "Prefix of job name: ", false, 10));
    	
    }
}

