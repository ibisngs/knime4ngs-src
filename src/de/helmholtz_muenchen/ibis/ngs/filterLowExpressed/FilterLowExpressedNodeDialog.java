package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

/**
 * <code>NodeDialog</code> for the "FilterLowExpressed" Node.
 * 
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelInteger SET_KEEP_READS	= new SettingsModelInteger(FilterLowExpressedNodeModel.CFGKEY_KEEP_READS, FilterLowExpressedNodeModel.DEFAULT_KEEP_READS);
    private final SettingsModelDouble SET_KEEP_FRACTION	= new SettingsModelDouble(FilterLowExpressedNodeModel.CFGKEY_KEEP_FRACTION, FilterLowExpressedNodeModel.DEFAULT_KEEP_FRATION);
    private final SettingsModelBoolean SET_BOTH_SEP		= new SettingsModelBoolean(FilterLowExpressedNodeModel.CFGKEY_BOTH_SEPERATE, FilterLowExpressedNodeModel.DEFAULT_BOTH_SEP);
 
	/**
	 * Constructor
	 */
    protected FilterLowExpressedNodeDialog() {
    	super();
    	// create the components
    	DialogComponentNumber keepReads = new DialogComponentNumber(SET_KEEP_READS, "minimum read number", 1);
    	DialogComponentNumber keepFraction = new DialogComponentNumber(SET_KEEP_FRACTION, "fraction of samples", 0.01);
    	DialogComponentBoolean bothSep = new DialogComponentBoolean(SET_BOTH_SEP, "test for both conditions separately");
    	
    	this.addDialogComponent(keepReads);
    	this.addDialogComponent(keepFraction);
    	this.addDialogComponent(bothSep);
    }
}

