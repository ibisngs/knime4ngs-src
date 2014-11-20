package de.helmholtz_muenchen.ibis.ngs.edgeR;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "EdgeR" Node.
 * 
 * @author Michael Kluge
 */
public class EdgeRNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelString SET_CORRECTION	= new SettingsModelString(EdgeRNodeModel.CFGKEY_CORRECTION_METHOD, EdgeRNodeModel.DEFAULT_CORRECTION_METHOD);
    private final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(EdgeRNodeModel.CFGKEY_NORMALIZE_METHOD_FACTOR, EdgeRNodeModel.DEFAULT_NORMALIZE_METHOD_FACTOR);

    /**
	 * Constructor
	 */
    protected EdgeRNodeDialog() {
    	super();
    	// create the components
    	DialogComponentStringSelection cor = new DialogComponentStringSelection(SET_CORRECTION, "p.value correction method", EdgeRNodeModel.CORRECTION_METHODS, false);
    	DialogComponentStringSelection normFac = new DialogComponentStringSelection(SET_NORM_FACTOR, "method for calculation of normalization factors", EdgeRNodeModel.NORM_FACTORS, false);
    	
    	this.addDialogComponent(normFac);
    	this.addDialogComponent(cor);
    }
}

