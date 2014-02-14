package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.utils.FastaFileNameFilter;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * <code>NodeDialog</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 * 
 * @author Michael Kluge
 */
public class FastaSelectorNodeDialog extends DefaultNodeSettingsPane {

	private static final int VISIBLE_ITEMS = 10;
	private static final String NO_SELECTION_MADE = "--- nothing selected yet ---";
	
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FASTA_DIR 				= FastaSelectorNodeModel.getSettingsModelString(FastaSelectorNodeModel.CFGKEY_FASTA_DIR);
    private final SettingsModelString SET_FASTA_FILE 				= FastaSelectorNodeModel.getSettingsModelString(FastaSelectorNodeModel.CFGKEY_FASTA_FILE);
    private final SettingsModelStringArray SET_FASTA_LIST_DISPLAY	= new SettingsModelStringArray(FastaSelectorNodeModel.CFGKEY_FASTA_LIST_DISPLAY, new String[0]);

	// components which must be accest inside a event handler
	private final DialogComponentButton DC_FASTA_DIR_BUTTON				= new DialogComponentButton("add all fasta files of this folder");
	private final DialogComponentButton DC_FASTA_FILE_BUTTON			= new DialogComponentButton("add selected fasta file");
	private final DialogComponentButton DC_FASTA_REMOVE_BUTTON			= new DialogComponentButton("remove selected fasta files");
	private final DialogComponentStringListSelection DC_FASTA_DISPLAY 	= new DialogComponentStringListSelection(SET_FASTA_LIST_DISPLAY, "files: ", NO_SELECTION_MADE);
    
	// storage 
	private final HashSet<String> FASTA_FILES = new HashSet<String>();
	
	 // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastaSelectorNodeDialog.class);
    
    @SuppressWarnings("deprecation")
	protected FastaSelectorNodeDialog() {
        super();
        
        // create dialogs
        DialogComponentFileChooser dcFastaDir 	= new DialogComponentFileChooser(SET_FASTA_DIR, "his_id_FASTA_DIR", 0, true);
        DialogComponentFileChooser dcFastaFile 	= new DialogComponentFileChooser(SET_FASTA_FILE, "his_id_FASTA_FILE", 0);

        // set title text
        dcFastaDir.setBorderTitle("path to fasta folder");
        dcFastaFile.setBorderTitle("path to fasta file");
    	
        // configure visible items 
        DC_FASTA_DISPLAY.setVisibleRowCount(VISIBLE_ITEMS);
        
        // disable both buttons at first time
        DC_FASTA_DIR_BUTTON.setEnabled(false);
        DC_FASTA_FILE_BUTTON.setEnabled(false);
        DC_FASTA_REMOVE_BUTTON.setEnabled(false);
               
        // add elements
        createNewGroup("add fasta files");
        addDialogComponent(dcFastaFile);
        addDialogComponent(DC_FASTA_FILE_BUTTON);
        addDialogComponent(dcFastaDir);
        addDialogComponent(DC_FASTA_DIR_BUTTON);
        
        createNewGroup("currently selected fasta files");
        addDialogComponent(DC_FASTA_DISPLAY);  
        addDialogComponent(DC_FASTA_REMOVE_BUTTON);
       
        // add action listener
        DC_FASTA_DIR_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addFastaFiles(SET_FASTA_DIR.getStringValue());
				SET_FASTA_DIR.setStringValue(FastaSelectorNodeModel.DEFAULT_FASTA_DIR);
				DC_FASTA_DIR_BUTTON.setEnabled(false);
			}
        });
        
        DC_FASTA_FILE_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addFastaFile(SET_FASTA_FILE.getStringValue());
				SET_FASTA_FILE.setStringValue(FastaSelectorNodeModel.DEFAULT_FASTA_FILE);
				DC_FASTA_FILE_BUTTON.setEnabled(false);
			}
        });
        
        DC_FASTA_REMOVE_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				removeFastaFile(SET_FASTA_LIST_DISPLAY.getStringArrayValue());
			}
        });
        
        // enable button on change of textbox
        SET_FASTA_DIR.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				if(SET_FASTA_DIR.getStringValue().length() > 0)
					DC_FASTA_DIR_BUTTON.setEnabled(true);
			}
        });
        SET_FASTA_FILE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				if(SET_FASTA_FILE.getStringValue().length() > 0)
					DC_FASTA_FILE_BUTTON.setEnabled(true);
			}
        });
    }
    
    @Override
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.FASTA_FILES.clear();
    	// check, if data is set
        if (settings.containsKey(FastaSelectorNodeModel.CFGKEY_FASTA_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(FastaSelectorNodeModel.CFGKEY_FASTA_LIST))
					this.FASTA_FILES.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }

    
    @Override
    public void saveAdditionalSettingsTo(NodeSettingsWO settings) {
    	// save the hash to the key FastaSelectorNodeModel.CFGKEY_FASTA_LIST
    	settings.addStringArray(FastaSelectorNodeModel.CFGKEY_FASTA_LIST, FASTA_FILES.toArray(new String[FASTA_FILES.size()]));
    }

    
    /**
     * adds all fasta files in the folder
     * @param dirname
     */
    private void addFastaFiles(String dirname) {
    	// get fasta files from folder
    	ArrayList<String> files = new ArrayList<String>();
    	files = IO.getFilesOfFolder(dirname, new FastaFileNameFilter());

    	// add the fasta files
    	for(String f : files)
    		addFastaFile(f);
    }
    
    /**
     * adds a single fasta file
     * @param file
     */
    @SuppressWarnings("deprecation")
	private void addFastaFile(String filename) {
    	// check, if file is there
    	File f = new File(filename);
    	if(!f.isFile() || !f.exists())
    		return;
    	
    	// check for file ending
    	if(!new FastaFileNameFilter().accept(null, filename))
    		return;
    	
    	// check, if fasta file is valid
    	if(!FileValidator.checkFastaFormat(filename))
    		return;
    	
    	// check, if filename is already stored
    	if(FASTA_FILES.contains(filename))
    		return;
    	
    	// add new value
    	FASTA_FILES.add(filename);
    	// update the list
    	DC_FASTA_DISPLAY.replaceListItems(FASTA_FILES, new String[0]);
    	DC_FASTA_DISPLAY.setVisibleRowCount(FASTA_FILES.size());
    	
    	// ensure that button is enabled
		DC_FASTA_REMOVE_BUTTON.setEnabled(true);
		DC_FASTA_DISPLAY.setVisibleRowCount(VISIBLE_ITEMS);
    }
    
    /**
     * removes filesnames from the selection
     * @param filenames
     */
    @SuppressWarnings("deprecation")
	private void removeFastaFile(String[] filenames) {
    	for(String f : filenames)
    		FASTA_FILES.remove(f);
    	
    	// update the list
    	// disable button if needed
		if(FASTA_FILES.size() == 0) {
			DC_FASTA_REMOVE_BUTTON.setEnabled(false);
			ArrayList<String> empty = new ArrayList<String>();
			empty.add(NO_SELECTION_MADE);
			DC_FASTA_DISPLAY.replaceListItems(empty, (String[]) null);
		}
		else
			DC_FASTA_DISPLAY.replaceListItems(FASTA_FILES, (String[]) null);
		
		// reset the number of rows which are displayed
		DC_FASTA_DISPLAY.setVisibleRowCount(VISIBLE_ITEMS);
	}
}
