package de.helmholtz_muenchen.ibis.ngs.vcfSampler;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VCFSampler" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class VCFSamplerNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelInteger ctrls = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CTRLS, 100);
	private final SettingsModelInteger buffer = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_BUFFER, 100);
	private final SettingsModelInteger cases = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CASES, 100);
	private final SettingsModelString def = new SettingsModelString(VCFSamplerNodeModel.CFGKEY_DEF,"10/2");
	
    /**
     * New pane for configuring the VCFSampler node.
     */
    protected VCFSamplerNodeDialog() {

    	addDialogComponent(new DialogComponentNumber(cases, "Number of cases", 100));
    	addDialogComponent(new DialogComponentNumber(ctrls, "Number of controls", 100));
    	addDialogComponent(new DialogComponentString(def, "Define variant frequency changes"));
    	addDialogComponent(new DialogComponentNumber(buffer, "Buffer size", 50));
    }
}

