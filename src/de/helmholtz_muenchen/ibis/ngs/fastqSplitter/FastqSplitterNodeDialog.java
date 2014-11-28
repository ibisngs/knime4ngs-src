package de.helmholtz_muenchen.ibis.ngs.fastqSplitter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumberEdit;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;
import org.knime.core.util.Pair;


/**
 * <code>NodeDialog</code> for the "FastqSplitter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Michael Kluge
 */
public class FastqSplitterNodeDialog extends DefaultNodeSettingsPane {
	
	private static final int VISIBLE_ITEMS = 10;
	private static final String NO_SELECTION_MADE = "--- no range added yet ---";
	private static final String ALLOWED_CHARS = "[A-Za-z0-9\\.\\-_]";

	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelInteger SET_MIN				= new SettingsModelInteger(FastqSplitterNodeModel.CFGKEY_MIN_LENGTH, FastqSplitterNodeModel.DEFAULT_MIN_LENGTH);
    private final SettingsModelInteger SET_MAX				= new SettingsModelInteger(FastqSplitterNodeModel.CFGKEY_MAX_LENGTH, FastqSplitterNodeModel.DEFAULT_MAX_LENGTH);
    private final SettingsModelString SET_NAME				= new SettingsModelString(FastqSplitterNodeModel.CFGKEY_BIN_NAME, FastqSplitterNodeModel.DEFAULT_BIN_NAME);
    private final SettingsModelStringArray SET_LIST_DISPLAY	= new SettingsModelStringArray(FastqSplitterNodeModel.CFGKEY_RANGES, new String[0]);
    private final SettingsModelString SET_OUTPUT_FOLDER		= new SettingsModelString(FastqSplitterNodeModel.CFGKEY_OUTPUT_FOLDER, FastqSplitterNodeModel.DEFAULT_OUTPUT_FOLDER);
    
	// components which must be accessible inside a event handler or somewhere else
	private final DialogComponentButton DC_ADD_BUTTON				= new DialogComponentButton("Add split range!");
	private final DialogComponentButton DC_REMOVE_BUTTON			= new DialogComponentButton("Remove selected split range!");
	private final DialogComponentButton DC_CLEAR_BUTTON				= new DialogComponentButton("Clear all set ranges!");
	private final DialogComponentStringListSelection DC_DISPLAY 	= new DialogComponentStringListSelection(SET_LIST_DISPLAY, "added ranges: ", NO_SELECTION_MADE);
	private final DialogComponentString DC_NAME 					= new DialogComponentString(SET_NAME, "name for the split range");
	private final DialogComponentNumberEdit DC_MIN					= new DialogComponentNumberEdit(SET_MIN, "minimum length");
	private final DialogComponentNumberEdit DC_MAX					= new DialogComponentNumberEdit(SET_MAX, "maximum length");
	private final DialogComponentFileChooser DC_OUTPUT			 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_split_outputFolder", 0, true);
	// storage 
	private final HashMap<String, Pair<Integer, Integer>> RANGES = new HashMap<String, Pair<Integer, Integer>>();
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastqSplitterNodeDialog.class);
	
    /**
     * New pane for configuring the FastqSplitter node.
     */
    protected FastqSplitterNodeDialog() {
    	
    	this.createNewGroup("add split range");
    	this.setHorizontalPlacement(true);
    	this.addDialogComponent(DC_NAME);
    	this.addDialogComponent(DC_MIN);
    	this.addDialogComponent(DC_MAX);
    	this.setHorizontalPlacement(false);
    	this.addDialogComponent(DC_ADD_BUTTON);

        // configure visible items 
        DC_DISPLAY.setVisibleRowCount(VISIBLE_ITEMS);
		ArrayList<String> empty = new ArrayList<String>();
		empty.add(NO_SELECTION_MADE);
		DC_DISPLAY.replaceListItems(empty, NO_SELECTION_MADE);

		this.createNewGroup("Ranges");
		this.addDialogComponent(DC_DISPLAY);
		this.addDialogComponent(DC_REMOVE_BUTTON);
		this.addDialogComponent(DC_CLEAR_BUTTON);
		
		this.createNewGroup("OutputSettings");
		this.addDialogComponent(DC_OUTPUT);
        
		DC_ADD_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addRange(SET_NAME.getStringValue(), SET_MIN.getIntValue(), SET_MAX.getIntValue());
			}
        });
        
        DC_REMOVE_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				removeRange(SET_LIST_DISPLAY.getStringArrayValue());
			}
        });
        
        DC_CLEAR_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				clearRange();
			}
        });
    }
    
    @Override
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.RANGES.clear();
    	// check, if data is set
        if (settings.containsKey(FastqSplitterNodeModel.CFGKEY_RANGES_NAMES)) {
        	try {
        		String[] names = settings.getStringArray(FastqSplitterNodeModel.CFGKEY_RANGES_NAMES);
        		int[] min = settings.getIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MIN);
        		int[] max = settings.getIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MAX);
        		// add the values
				for(int i = 0; i < names.length; i++) {
					this.addRange(names[i], min[i], max[i]);
				}
				
				// update the list
		    	this.DC_DISPLAY.replaceListItems(getNamesForList(), (String[]) null);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }


    @Override
    public void saveAdditionalSettingsTo(NodeSettingsWO settings) {
    	ArrayList<String> n = new ArrayList<String>();
    	int min[] = new int[this.RANGES.size()];
    	int max[] = new int[this.RANGES.size()];

    	int i = 0;
    	for(String name : this.RANGES.keySet()) {
    		Pair<Integer, Integer> p = this.RANGES.get(name);
    		n.add(name);
    		min[i] = p.getFirst();
    		max[i] = p.getSecond();
    		i++;
    	}
    	// save the hash to the keys
    	settings.addStringArray(FastqSplitterNodeModel.CFGKEY_RANGES_NAMES, n.toArray(new String[n.size()]));
    	settings.addIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MIN, min);
    	settings.addIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MAX, max);
    }
    
     
    /**
     * adds a single Range
     * @param rangename
     * @param min
     * @param max
     */
	private void addRange(String rangename, int min, int max) {
      	// check, if rangename is already stored
    	if(this.RANGES.containsKey(rangename))
    		return;
    	
    	// check, if range is valid
    	if(!isRangeValid(rangename, min, max))
    		return;

    	// add new value
    	this.RANGES.put(rangename, new Pair<Integer, Integer>(min, max));
    	// update the list
    	this.DC_DISPLAY.replaceListItems(getNamesForList(), (String[]) null);

    	// ensure that button is enabled
    	//this.DC_REMOVE_BUTTON.getModel().setEnabled(true);
    	//this.DC_CLEAR_BUTTON.getModel().setEnabled(true);
    }
    
    /**
     * removes rangenames from the selection
     * @param rangenames
     */
	private void removeRange(String[] rangenames) {
    	for(String n : rangenames)
    		this.RANGES.remove(n.split(" \\[")[0]);
    	
    	// update the list
    	// disable button if needed
		if(this.RANGES.size() == 0) {
			//this.DC_REMOVE_BUTTON.getModel().setEnabled(false);
			//this.DC_CLEAR_BUTTON.getModel().setEnabled(false);
			ArrayList<String> empty = new ArrayList<String>();
			empty.add(NO_SELECTION_MADE);
			this.DC_DISPLAY.replaceListItems(empty, NO_SELECTION_MADE);
		}
		else
			this.DC_DISPLAY.replaceListItems(this.getNamesForList(), (String[]) null);
	}
	
	/**
	 * clears all of the stuff
	 */
	private void clearRange() {
		this.RANGES.clear();
		//this.DC_REMOVE_BUTTON.getModel().setEnabled(false);
		//this.DC_CLEAR_BUTTON.getModel().setEnabled(false);
		
		ArrayList<String> empty = new ArrayList<String>();
		empty.add(NO_SELECTION_MADE);
		
		this.DC_DISPLAY.replaceListItems(empty, NO_SELECTION_MADE);
	}
	
             
    /**
     * checks, if a Range is valid
     * @param rangename
     * @param min
     * @param max
     */
    public boolean isRangeValid(String rangename, int min, int max) {
    	if(min >= max) {
    		LOGGER.error("Min must be greater than max value.");
    		return false;
    	}
    	
    	if(rangename.replaceAll(ALLOWED_CHARS, "").length() > 0) {
    		LOGGER.error("File name can only contain: " + ALLOWED_CHARS);
    		return false;
    	}
    	return true;
    }
    
    /**
     * returns the values for the list
     * @return
     */
    private ArrayList<String> getNamesForList() {
    	ArrayList<String> n = new ArrayList<String>();
    	for(String name : this.RANGES.keySet()) {
    		Pair<Integer, Integer> p = this.RANGES.get(name);
    		n.add(name + " [" + p.getFirst() + ", " + p.getSecond() + "]");
    	}
    	return n;
    }
}
