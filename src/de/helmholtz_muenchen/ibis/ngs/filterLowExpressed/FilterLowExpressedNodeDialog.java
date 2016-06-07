package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "FilterLowExpressed" Node.
 * 
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelInteger SET_KEEP_READS	= new SettingsModelInteger(FilterLowExpressedNodeModel.CFGKEY_KEEP_READS, FilterLowExpressedNodeModel.DEFAULT_KEEP_READS);
    private final SettingsModelDouble SET_KEEP_FRACTION	= new SettingsModelDouble(FilterLowExpressedNodeModel.CFGKEY_KEEP_FRACTION, FilterLowExpressedNodeModel.DEFAULT_KEEP_FRATION);
    private final SettingsModelBoolean SET_BOTH_SEP		= new SettingsModelBoolean(FilterLowExpressedNodeModel.CFGKEY_BOTH_SEPERATE, FilterLowExpressedNodeModel.DEFAULT_BOTH_SEP);
    private final SettingsModelString SET_MODE 			= new SettingsModelString(FilterLowExpressedNodeModel.CFGKEY_MODE, FilterLowExpressedNodeModel.DEFAULT_MODE);
    
	/**
	 * Constructor
	 */
    protected FilterLowExpressedNodeDialog() {
    	super();
    	// create the components
    	DialogComponentNumber keepReads = new DialogComponentNumber(SET_KEEP_READS, "Minimum read number", 1);
    	DialogComponentNumber keepFraction = new DialogComponentNumber(SET_KEEP_FRACTION, "Fraction of samples", 0.01);
    	DialogComponentBoolean bothSep = new DialogComponentBoolean(SET_BOTH_SEP, "Filter both conditions separately");

		DialogComponentButtonGroup mode = new DialogComponentButtonGroup(SET_MODE, false, "Filtering Mode", FilterLowExpressedNodeModel.DEFAULT_MODE, "Per sample");

    	this.addDialogComponent(mode);
    	this.addDialogComponent(keepReads);
    	this.addDialogComponent(keepFraction);
    	this.addDialogComponent(bothSep);
    	
		// check for changes
    	this.SET_MODE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(SET_MODE.getStringValue().equals(FilterLowExpressedNodeModel.DEFAULT_MODE)) {
					SET_KEEP_FRACTION.setEnabled(false);
				}
				else {
					SET_KEEP_FRACTION.setEnabled(true);
				}
			}
        });
    }
}

