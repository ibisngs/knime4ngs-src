package de.helmholtz_muenchen.ibis.ngs.DESeq;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "DESeq" Node.
 * 
 * @author Michael Kluge
 */
public class DESeqNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelString SET_METHOD	= new SettingsModelString(DESeqNodeModel.CFGKEY_METHOD, DESeqNodeModel.DEFAULT_METHOD);
    private final SettingsModelString SET_SHEARING	= new SettingsModelString(DESeqNodeModel.CFGKEY_SHEARING, DESeqNodeModel.DEFAULT_SHEARING);
    @SuppressWarnings("unused")
	private final SettingsModelString SET_VS			= new SettingsModelString(DESeqNodeModel.CFGKEY_VS, DESeqNodeModel.DEFAULT_VS);

	/**
	 * Constructor
	 */
    protected DESeqNodeDialog() {
    	super();
    	// create the components
    	DialogComponentStringSelection method = new DialogComponentStringSelection(SET_METHOD, "empirical dispersion calculation", DESeqNodeModel.METHODS, false);
    	DialogComponentStringSelection shearing = new DialogComponentStringSelection(SET_SHEARING, "sharing mode", DESeqNodeModel.SHEARING, false);
    	
    	this.addDialogComponent(method);
    	this.addDialogComponent(shearing);
    }
}

