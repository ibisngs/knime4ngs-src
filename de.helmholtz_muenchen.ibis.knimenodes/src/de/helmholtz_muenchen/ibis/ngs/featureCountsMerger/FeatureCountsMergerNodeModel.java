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
package de.helmholtz_muenchen.ibis.ngs.featureCountsMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.IntCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

/**
 * This is the model implementation of FeatureCountsMerger.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsMergerNodeModel extends HTExecutorNodeModel {
    
	public static final String ID = "ID";
	public static final String TAB = "\t";
	public static final String ENDING_SPLIT = "\\.";
	
	// keys for SettingsModels
	protected static final String CFGKEY_OUTPUT_FILE 	= "OutputFile";
	protected static final String CFGKEY_REMOVE_ENDING 	= "RemoveEnding";
	protected static final String CFGKEY_REMOVE_PATH 	= "RemoevPath";
	
    // initial default values for SettingsModels    
	protected static final String DEFAULT_OUTPUT_FILE 	= "";
	protected static final boolean DEFAULT_REMOVE_ENDING= true;
	protected static final boolean DEFAULT_REMOVE_PATH 	= true;
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(CFGKEY_OUTPUT_FILE, DEFAULT_OUTPUT_FILE);
    private final SettingsModelBoolean SET_REMOVE_ENDING = new SettingsModelBoolean(CFGKEY_REMOVE_ENDING, DEFAULT_REMOVE_ENDING);
    private final SettingsModelBoolean SET_REMOVE_PATH = new SettingsModelBoolean(CFGKEY_REMOVE_PATH, DEFAULT_REMOVE_PATH);

    /**
     * Constructor for the node model.
     */
    protected FeatureCountsMergerNodeModel() {
        super(1, 1);
        addSetting(SET_OUTPUT_FILE);
    	addSetting(SET_REMOVE_ENDING);
    	addSetting(SET_REMOVE_PATH);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    //overrides HTE execution method as no external tool is used
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	
    	String outfile = IO.processFilePath(this.SET_OUTPUT_FILE.getStringValue());
    	if(CompatibilityChecker.inputFileNotOk(outfile, false)) {
    		String in = inData[0].iterator().next().getCell(0).toString();
    		outfile = new File(in).getParent() + File.separator + "merged.featureCounts";
    	}
    	
		boolean do_overwrite;
		if(m_use_pref.getBooleanValue()) {
			do_overwrite = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.OVERWRITE);
		} else {
			do_overwrite = m_overwrite.getBooleanValue();
		}
		
		//abort execution if node shall not overwrite existing outfiles
		if(!do_overwrite && Files.exists(Paths.get(outfile))) {
			throw new UnsuccessfulExecutionException("Execution aborted as outfile "+outfile+" exists yet! Rename/move/delete existing outfile or allow overwriting of existing outfiles.");
		}
    	
    	// get featureCount files
    	ArrayList<String> dataFiles = new ArrayList<String>();
    	for(CloseableRowIterator it = inData[0].iterator(); it.hasNext(); ) {
    		dataFiles.add(it.next().getCell(0).toString());
    	}
    	
    	// data gets stored here
    	HashMap<String, HashMap<String, String>> counts = new HashMap<String, HashMap<String, String>>();
    	ArrayList<String> sampleNames = new ArrayList<String>();
    	
    	// run though all files
    	for(String fname : dataFiles) {
    		File f = new File(fname);
    		// test if correct file is given
    		if(!f.isFile() || !f.canRead())
    			throw new IllegalArgumentException("File '" + fname + "' was not found.");
    		
    		String sampleName = f.getAbsolutePath();
    		// remove path, if variable is set
    		if(this.SET_REMOVE_PATH.getBooleanValue()) 
    			sampleName = f.getName();
    		// remove ending (stuff after last ".") if variable is set
    		if(this.SET_REMOVE_ENDING.getBooleanValue()) {
    			String tmp[] = sampleName.split(ENDING_SPLIT);
    			sampleName = tmp[0];
//    			for(int i = 1; i < tmp.length-1; i++) {
//    				sampleName = sampleName + tmp[i];
//    			}
    		}

    		sampleNames.add(sampleName);
    		
    		// open file
    		try {
    			BufferedReader r = new BufferedReader(new FileReader(f));
    			String line = null;
    			String name = null;
    			String count = null;
    			boolean headerLine = true;
    			
    			// read lines
    			while((line = r.readLine()) != null) {
    				// try to split line
    				String tmp[] = line.split(TAB);
    				if(tmp.length == 7) {
    					// replace white spaces at start and end, if some are there
    					name = tmp[0].trim();
    					count = tmp[6].trim();		
    					
    					// store data (skip first line)
    					if(!headerLine) {
    						// check if HashMap for gene is there
    						if(!counts.containsKey(name))
    							counts.put(name, new HashMap<String, String>());
    						
    						// get hashmap and store data
    						counts.get(name).put(sampleName, count);
    					}
    					headerLine = false;
    				}	
    			}
    			// close the file
    			r.close();
    		}
    		catch(Exception e) {
    			e.printStackTrace(); 
    		}
    	}

    	// prepare output file and BufferedDataTable
    	DataTableSpec table = getDataOutSpec(sampleNames);
        BufferedDataContainer cont = exec.createDataContainer(table);
    	
//    	System.out.println(outfile);
    	
//    	new File(outfile).delete();
    	BufferedWriter outfileBW = new BufferedWriter(new FileWriter(outfile));
    	
    	// write header
    	String sn;
    	outfileBW.write(ID);
    	outfileBW.write(TAB);
    	for(Iterator<String> it = sampleNames.iterator(); it.hasNext(); ) {
    		outfileBW.write(it.next());
    		if(it.hasNext()) 
    			outfileBW.write(TAB);
    	}
    	outfileBW.newLine();
    	outfileBW.flush();
    	
    	// run though all files and get counts
    	for(Iterator<String> itNames = counts.keySet().iterator(); itNames.hasNext(); ) {
    		String name = itNames.next();
    		HashMap<String, String> c = counts.get(name);
    		
    		// write name of feature
    		outfileBW.write(name);
        	outfileBW.write(TAB);
        	
        	// write counts
        	DataCell[] cellArray = new DataCell[sampleNames.size()];
        	int i = 0;
    		for(Iterator<String> it = sampleNames.iterator(); it.hasNext(); ) {
    			int v = 0;
    			sn = it.next();
    			if(c.containsKey(sn)) {
    				outfileBW.write(c.get(sn));
    				v = Integer.parseInt(c.get(sn));
    			}
    			else
    				outfileBW.write("0");
    					
    			if(it.hasNext()) 
        			outfileBW.write(TAB);
    			
    	    	// write output to table
    	    	IntCell cell = new IntCell(v);
    	    	cellArray[i] = cell;
    	    	i++;

    		}
    		// add row to table
	    	DefaultRow r = new DefaultRow(name, cellArray);
	    	cont.addRowToTable(r);
	    	
    		// make new line
    		if(itNames.hasNext())
    			outfileBW.newLine();
    	}
    	
    	outfileBW.flush();
    	// close outfile
    	outfileBW.close();

    	// finalize table
    	cont.close();
    	BufferedDataTable out = cont.getTable();
        return new BufferedDataTable[]{ out };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return new DataTableSpec[]{ null };
    }
    
    /**
     * get specs of table
     * @param sampleNames
     * @return
     */
    private DataTableSpec getDataOutSpec(ArrayList<String> sampleNames) {
    	DataType[] types = new DataType[sampleNames.size()];
    	// add types for int
    	for(int i = 0; i < sampleNames.size(); i++)
    		types[i] = DataType.getType(IntCell.class);
    	
    	DataColumnSpec[] specs = DataTableSpec.createColumnSpecs(sampleNames.toArray(new String[0]), types);
    	return(new DataTableSpec(specs));
	}
}

