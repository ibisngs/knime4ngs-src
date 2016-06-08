package de.helmholtz_muenchen.ibis.ngs.limma;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Limma" Node.
 * 
 * @author Michael Kluge
 */
public class LimmaNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelString SET_CORRECTION	= new SettingsModelString(LimmaNodeModel.CFGKEY_CORRECTION_METHOD, LimmaNodeModel.DEFAULT_CORRECTION_METHOD);
    private final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(LimmaNodeModel.CFGKEY_NORMALIZE_METHOD_FACTOR, LimmaNodeModel.DEFAULT_NORMALIZE_METHOD_FACTOR);
    private final SettingsModelString SET_METHOD_CPM	= new SettingsModelString(LimmaNodeModel.CFGKEY_NORMALIZE_METHOD_CPM, LimmaNodeModel.DEFAULT_NORMALIZE_METHOD_CPM);
    
	/**
	 * Constructor
	 */
    protected LimmaNodeDialog() {
    	super();
    	// create the components
    	DialogComponentStringSelection cor = new DialogComponentStringSelection(SET_CORRECTION, "P-value correction method", LimmaNodeModel.CORRECTION_METHODS, false);
    	DialogComponentStringSelection normFac = new DialogComponentStringSelection(SET_NORM_FACTOR, "Method for calculation of normalization factors", LimmaNodeModel.NORM_FACTORS, false);
    	DialogComponentStringSelection normCPM = new DialogComponentStringSelection(SET_METHOD_CPM, "Method for CPM normalization", LimmaNodeModel.NORM_CPM, false);
    	
    	this.addDialogComponent(normFac);
    	this.addDialogComponent(normCPM);
    	this.addDialogComponent(cor);
    }
}

