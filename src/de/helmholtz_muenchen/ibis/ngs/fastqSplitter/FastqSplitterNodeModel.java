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

package de.helmholtz_muenchen.ibis.ngs.fastqSplitter;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.util.Pair;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * This is the model implementation of FastqSplitter.
 * 
 *
 * @author Michael Kluge
 */
public class FastqSplitterNodeModel extends SettingsStorageNodeModel {
	
	// TODO: add node description
	
	// keys for SettingsModels
    protected static final String CFGKEY_MIN_LENGTH = "SplitterMinLength";
    protected static final String CFGKEY_MAX_LENGTH = "SplitterMaxLength";
    protected static final String CFGKEY_BIN_NAME	= "SplitterBinName";
    protected static final String CFGKEY_RANGES = "SplitterListBins";
    protected static final String CFGKEY_RANGES_NAMES = "SplitterListNames";
    protected static final String CFGKEY_RANGES_MIN = "SplitterListMin";
    protected static final String CFGKEY_RANGES_MAX = "SplitterListMax";
    protected static final String CFGKEY_OUTPUT_FOLDER = "SplitterOutputFolder";
    
    // initial default values for SettingsModels
    protected static final int DEFAULT_MIN_LENGTH = 0;
    protected static final int DEFAULT_MAX_LENGTH = 100;
    protected static final String DEFAULT_BIN_NAME = "NameOfSplit";
    protected static final String DEFAULT_OUTPUT_FOLDER = "";
    
    public static final String OUTPUT_NAME_RANGENAME = "RangeName";
    public static final String OUTPUT_NAME_MIN = "RangeMinLength";
    public static final String OUTPUT_NAME_MAX = "RangeMaxLength";
    public static final String OUTPUT_NAME_PATH = "Path2Fastq";
    public static final String ENDING = ".fastq";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelInteger SET_MIN	= new SettingsModelInteger(CFGKEY_MIN_LENGTH, DEFAULT_MIN_LENGTH);
    private final SettingsModelInteger SET_MAX	= new SettingsModelInteger(CFGKEY_MAX_LENGTH, DEFAULT_MAX_LENGTH);
    private final SettingsModelString SET_NAME	= new SettingsModelString(CFGKEY_BIN_NAME, DEFAULT_BIN_NAME);
    private final SettingsModelString SET_OUTPUT_FOLDER	= new SettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
    
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastqSplitterNodeModel.class);
	private boolean hasConfigureOpendOnce 	= false; // true, if configure was opend once
	private final HashMap<String, Pair<Integer, Integer>> RANGES = new HashMap<String, Pair<Integer, Integer>>();
	
	
    /**
     * Constructor for the node model.
     */
    protected FastqSplitterNodeModel() {
        super(1, 1);
        this.addSetting(SET_MIN);
		this.addSetting(SET_MAX);
		this.addSetting(SET_NAME);
		this.addSetting(SET_OUTPUT_FOLDER);
    }
	
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// check, if some files where selected
    	if(hasConfigureOpendOnce && this.RANGES.size() == 0)
    		throw new InvalidSettingsException("Add at least one range.");
    	
    	// check, for valid output folder
    	try {
    		new File(this.SET_OUTPUT_FOLDER.getStringValue()).getCanonicalPath();
    	}
    	catch(Exception e) {
    		throw new InvalidSettingsException("Enter a valid output folder path!");
    	}
    	
    	return new DataTableSpec[]{getDataOutSpec()};
    }
    
    /**
     * returns the output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec() {
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator(OUTPUT_NAME_RANGENAME, StringCell.TYPE);
    	DataColumnSpecCreator col2 = new DataColumnSpecCreator(OUTPUT_NAME_MIN, IntCell.TYPE);
    	DataColumnSpecCreator col3 = new DataColumnSpecCreator(OUTPUT_NAME_MAX, IntCell.TYPE);
    	DataColumnSpecCreator col4 = new DataColumnSpecCreator(OUTPUT_NAME_PATH, StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(), col2.createSpec(), col3.createSpec(), col4.createSpec()};
        
        return new DataTableSpec(cols);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	exec.setProgress(0.01);
    	
    	/************************* CHECK IF INPUT IS VALID ********/
    	String fastqFile = inData[0].iterator().next().getCell(0).toString();
    	File f = new File(fastqFile);
    	
    	if(!FileValidator.checkFastqFormat(f.getAbsolutePath()))
    		throw new InvalidSettingsException("'" + f.getAbsolutePath() + "' is not a valid FASTQ file.");
    	
    	/******************* PREPARE OUTPUT **********************/
    	DataTableSpec table = getDataOutSpec();
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	
    	// write output to table
    	File outf = new File(this.SET_OUTPUT_FOLDER.getStringValue());
    	int ii = 0;
    	for(String range : this.RANGES.keySet()) {
    		String fname = outf.getAbsolutePath() + File.separator + range + ENDING;
    		Pair<Integer, Integer> p = this.RANGES.get(range);
    		StringCell cl1 = new StringCell(range);
    		IntCell cl2 = new IntCell(p.getFirst());
    		IntCell cl3 = new IntCell(p.getSecond());
    		StringCell cl4 = new StringCell(new File(fname).getAbsolutePath());
	    	DataCell[] c = new DataCell[]{cl1, cl2, cl3, cl4};
	    	DefaultRow r = new DefaultRow("Row"+ii, c);
	    	cont.addRowToTable(r);
	    	ii++;
    	}
    	
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	exec.setProgress(0.02);
    	
    	/******************* PREPARE FILEHANDLES **********************/
    	HashSet<BufferedOutputStream> outputBuffers = new HashSet<BufferedOutputStream>();
    	HashMap<Integer, ArrayList<BufferedOutputStream>> length = new HashMap<Integer, ArrayList<BufferedOutputStream>>();
    	// create output folder
    	if(outf.isDirectory() || outf.mkdirs()) {
    		// create the filehandles
    		for(String range : this.RANGES.keySet()) {
    			String fname = outf.getAbsolutePath() + File.separator + range + ENDING;
    			BufferedOutputStream b = new BufferedOutputStream(new FileOutputStream(new File(fname)));
    			outputBuffers.add(b);
    			
    			// add this buffer to the complete range
    			Pair<Integer, Integer> p = this.RANGES.get(range);
    			for(Integer i = p.getFirst(); i <= p.getSecond(); i++) {
    				// check, if there is already an arrayList
    				if(!length.containsKey(i)) {
    					length.put(i, new ArrayList<BufferedOutputStream>());
    				}
    				// add entry
    				length.get(i).add(b);
    			}
    			
    		}
    		exec.setProgress(0.05);
    		
    		
    		/******************* OPEN FILE AND DOOOOOOO IT **********************/
    		BufferedReader br = null;
    		try {
    			String line, id, seq, desc, qual;
    			br = new BufferedReader(new FileReader(f.getAbsolutePath()));
    			int l = 0;
    			byte[] lineSep = System.lineSeparator().getBytes();
    			
    			while((line = br.readLine())!= null) {
    				if(line.length() > 0 && !line.startsWith("#")) { //if not an empty line or comment
    					if(line.startsWith("@")) {
    						id = line;
    						seq = br.readLine();
    						desc = br.readLine();
    						qual = br.readLine();
    						
    						// add it to the files
    						l = seq.length();
    						if(length.containsKey(l) && qual != null) {
    							for(BufferedOutputStream b : length.get(l)) {
    								// write the stuff
    								b.write(id.getBytes());
    								b.write(lineSep);
    								b.write(seq.getBytes());
    								b.write(lineSep);
    								b.write(desc.getBytes());
    								b.write(lineSep);
    								b.write(qual.getBytes());
    								b.write(lineSep);
    							}
    						}
    					}
    				}
    			}		
    		} catch(IOException e) {
    			e.printStackTrace();
    		} finally {
    			try {
    				if(br != null)
    					br.close();
    			} catch(IOException ex) {
    				ex.printStackTrace();
    			}
    		}
    		
    		/******************* !!!!!!!!!!!!!!!!!!!!!!!!! **********************/
    		// close all the output files
    		exec.setProgress(0.95);
    		for(BufferedOutputStream b : outputBuffers) {
    			b.flush();
    			b.close();
    		}

        	exec.setProgress(1.0); // we are ready! ;)
            return new BufferedDataTable[]{out};
    	}
    	else {
    		throw new InvalidSettingsException("Output folder '"+ outf.getAbsolutePath() +"' could not be created.");
    	}
    }
	
	
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	// clean the old data
    	this.RANGES.clear();
    	// check, if data is set
        if (settings.containsKey(FastqSplitterNodeModel.CFGKEY_RANGES_NAMES)) {
        	try {
        		String[] names = settings.getStringArray(FastqSplitterNodeModel.CFGKEY_RANGES_NAMES);
        		int[] min = settings.getIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MIN);
        		int[] max = settings.getIntArray(FastqSplitterNodeModel.CFGKEY_RANGES_MAX);
        		// add the values
				for(int i = 0; i < names.length; i++)
					this.RANGES.put(names[i], new Pair<Integer, Integer>(min[i], max[i]));
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }
    
    @Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);

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

	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		super.validateSettings(settings);
        // configure must have been opened or we won't be here
        hasConfigureOpendOnce = true;
	}
	
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {}

    @Override
    protected void reset() {
    	this.RANGES.clear();
    }
}

