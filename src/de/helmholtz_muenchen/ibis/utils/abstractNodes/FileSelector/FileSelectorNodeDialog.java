/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.utils.IO;

/**
 * <code>NodeDialog</code> for the abstract "FileSelector" Node.
 * This Node can be used to select multiple files.
 * 
 * @author Michael Kluge
 */
public abstract class FileSelectorNodeDialog extends DefaultNodeSettingsPane {

	private static final int VISIBLE_ITEMS = 10;
	private static final String NO_SELECTION_MADE = "--- nothing selected yet ---";
	
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FILE_DIR 					= new SettingsModelString(FileSelectorNodeModel.CFGKEY_FILE_DIR, FileSelectorNodeModel.DEFAULT_FILE_DIR);
    private final SettingsModelString SET_FILE_FILE 				= new SettingsModelString(FileSelectorNodeModel.CFGKEY_FILE_FILE, FileSelectorNodeModel.DEFAULT_FILE_FILE);
    private final SettingsModelString SET_NAME_REGEX				= new SettingsModelString(FileSelectorNodeModel.CFGKEY_REGEX, FileSelectorNodeModel.DEFAULT_REGEX);
    private final SettingsModelStringArray SET_FILE_LIST_DISPLAY	= new SettingsModelStringArray(FileSelectorNodeModel.CFGKEY_FILE_LIST_DISPLAY, new String[0]);
    
	// components which must be accessible inside a event handler or somewhere else
	private final DialogComponentButton DC_FILE_DIR_BUTTON				= new DialogComponentButton("add all " + getFiletypeName() + " files of this folder");
	private final DialogComponentButton DC_FILE_FILE_BUTTON				= new DialogComponentButton("add selected " + getFiletypeName() + " file");
	private final DialogComponentButton DC_FILE_REMOVE_BUTTON			= new DialogComponentButton("remove selected " + getFiletypeName() + " files");
	private final DialogComponentStringListSelection DC_FILE_DISPLAY 	= new DialogComponentStringListSelection(SET_FILE_LIST_DISPLAY, "files: ", NO_SELECTION_MADE);
    private final DialogComponentString DC_REGEX 						= new DialogComponentString(SET_NAME_REGEX, "filename regex filter");
    
	// storage 
	private final HashSet<String> FILE_FILES = new HashSet<String>();
	
	 // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FileSelectorNodeDialog.class);
    
	protected FileSelectorNodeDialog() {
        super();
        
        // create dialogs
        DialogComponentFileChooser dcFastaDir 	= new DialogComponentFileChooser(SET_FILE_DIR, "his_id_FILE_DIR", 0, true);
        DialogComponentFileChooser dcFastaFile 	= new DialogComponentFileChooser(SET_FILE_FILE, "his_id_FILE_FILE", 0);

        // set title text
        dcFastaDir.setBorderTitle("path to " + getFiletypeName() + " folder");
        dcFastaFile.setBorderTitle("path to " + getFiletypeName() + " file");
    	
        // configure visible items 
        DC_FILE_DISPLAY.setVisibleRowCount(VISIBLE_ITEMS);
		ArrayList<String> empty = new ArrayList<String>();
		empty.add(NO_SELECTION_MADE);
		DC_FILE_DISPLAY.replaceListItems(empty, NO_SELECTION_MADE);
		
		// TODO: get buttons working!!!! --> java.lang.AssertionError
        // disable both buttons at first time
        //DC_FILE_DIR_BUTTON.getModel().setEnabled(false);
        //DC_FILE_FILE_BUTTON.getModel().setEnabled(false);
        //DC_FILE_REMOVE_BUTTON.getModel().setEnabled(false);
               
        // add elements
        createNewGroup("add " + getFiletypeName() + " files");
        addDialogComponent(dcFastaFile);
        addDialogComponent(DC_FILE_FILE_BUTTON);
        
        addDialogComponent(dcFastaDir);
        setHorizontalPlacement(true);
        addDialogComponent(DC_REGEX);
        addDialogComponent(DC_FILE_DIR_BUTTON);
        setHorizontalPlacement(false);
        
        createNewGroup("currently selected " + getFiletypeName() + " files");
        addDialogComponent(DC_FILE_DISPLAY);  
        addDialogComponent(DC_FILE_REMOVE_BUTTON);
       
        // add action listener        
        DC_FILE_DIR_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addFiles(SET_FILE_DIR.getStringValue());
				SET_FILE_DIR.setStringValue(FileSelectorNodeModel.DEFAULT_FILE_DIR);
	//			DC_FILE_DIR_BUTTON.getModel().setEnabled(false);
			}
        });
        
        DC_FILE_FILE_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addFile(SET_FILE_FILE.getStringValue());
				SET_FILE_FILE.setStringValue(FileSelectorNodeModel.DEFAULT_FILE_FILE);
	//			DC_FILE_FILE_BUTTON.getModel().setEnabled(false);
			}
        });
        
        DC_FILE_REMOVE_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				removeFile(SET_FILE_LIST_DISPLAY.getStringArrayValue());
			}
        });
        
        // enable button on change of textbox
        SET_FILE_DIR.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
		//		if(SET_FILE_DIR.getStringValue().length() > 0)
		//			DC_FILE_DIR_BUTTON.getModel().setEnabled(true);
			}
        });
        SET_FILE_FILE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
		//		if(SET_FILE_FILE.getStringValue().length() > 0)
		//			DC_FILE_FILE_BUTTON.getModel().setEnabled(true);
			}
        });
    }
    
    @Override
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.FILE_FILES.clear();
    	// check, if data is set
        if (settings.containsKey(FileSelectorNodeModel.CFGKEY_FILE_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(FileSelectorNodeModel.CFGKEY_FILE_LIST))
					this.addFile(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }

    
    @Override
    public void saveAdditionalSettingsTo(NodeSettingsWO settings) {
    	// save the hash to the key FastaSelectorNodeModel.CFGKEY_FILE_LIST
    	settings.addStringArray(FileSelectorNodeModel.CFGKEY_FILE_LIST, FILE_FILES.toArray(new String[FILE_FILES.size()]));
    }
    
    /**
     * adds all files in the folder
     * @param dirname
     */
    private void addFiles(String dirname) {
    	// get files from folder
    	ArrayList<String> files = new ArrayList<String>();
    	files = IO.getFilesOfFolder(dirname, getFilenameFilter());
    	
    	// try to compile the REGEX
    	String regex = SET_NAME_REGEX.getStringValue();
    	Pattern p = null;
    	Matcher m = null;
    	try {
    		p = Pattern.compile(regex);
    	}
    	catch(Exception e) {
    		LOGGER.error("Pattern does not compile: " + e.getMessage());
    	}
    	
    	// add the files
    	for(String f : files) {
    		// pattern was not valid --> add all
    		if(p == null) 
    			addFile(f);
    		else {
    			m = p.matcher(new File(f).getName().replaceFirst(getFilenameEndRegex(), "")); //compile the pattern
    			// check, if pattern matches
    			if(m.matches())
    				addFile(f);
    		}
    	}
    }
    
    /**
     * adds a single file
     * @param file
     */
	private void addFile(String filename) {
    	// check, if file is there
    	File f = new File(filename);
    	if(!f.isFile() || !f.exists())
    		return;
    	
    	// check, if filename is already stored
    	if(FILE_FILES.contains(filename))
    		return;
    	
    	// check, if filename is valid
    	if(!isFileValid(f))
    		return;

    	// add new value
    	FILE_FILES.add(filename);
    	// update the list
    	DC_FILE_DISPLAY.replaceListItems(FILE_FILES, new String[0]);
    	
    	// ensure that button is enabled
		//DC_FILE_REMOVE_BUTTON.getModel().setEnabled(true);
    }
    
    /**
     * removes filenames from the selection
     * @param filenames
     */
	private void removeFile(String[] filenames) {
    	for(String f : filenames)
    		FILE_FILES.remove(f);
    	
    	// update the list
    	// disable button if needed
		if(FILE_FILES.size() == 0) {
		//	DC_FILE_REMOVE_BUTTON.getModel().setEnabled(false);
			ArrayList<String> empty = new ArrayList<String>();
			empty.add(NO_SELECTION_MADE);
			DC_FILE_DISPLAY.replaceListItems(empty, NO_SELECTION_MADE);
		}
		else
			DC_FILE_DISPLAY.replaceListItems(FILE_FILES, (String[]) null);
	}
	
    /**
     * must return a filename filter which is used when the folder option is used
     * can be overridden to extend to the function of just checking for the ending.
     * @return
     */
    public FilenameFilter getFilenameFilter() {
    	return new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				return name.matches(".*" + getFilenameEndRegex());
			}
    	};
    }
    
    
    /************************************** ABSTRACT METHODS **********************************************/
	/******************************************************************************************************/
    
    /**
     * Returns the filetype how it is named in the GUI
     * @return
     */
    protected abstract String getFiletypeName();
        
    /**
     * checks, if a file is valid for a specific filter
     * @param filename
     * @return
     */
    public abstract boolean isFileValid(File file);
    
    
    /**
     * returns a regex which matches the ends of the files
     * @return
     */
    public abstract String getFilenameEndRegex();
}
