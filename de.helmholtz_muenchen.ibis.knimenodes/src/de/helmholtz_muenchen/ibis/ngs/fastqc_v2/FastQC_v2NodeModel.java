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
package de.helmholtz_muenchen.ibis.ngs.fastqc_v2;

import java.io.File;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


/**
 * This is the model implementation of FastQC_v2.
 * 
 *
 * @author Paul Hager
 */
public class FastQC_v2NodeModel extends HTExecutorNodeModel {
	
	// the logger instance
    protected static final NodeLogger LOGGER = NodeLogger.getLogger(FastQC_v2NodeModel.class);
    
    // Node specific params
    public static final String CFGKEY_FASTQC = "Path2FastQC";
    public static final String CFGKEY_OUTFOLDER = "Path2OutFolder";
	public static final String CFGKEY_THREADS = "NumberThreads";
	public static final String CFGKEY_ADDITIONAL_OPTIONS = "AdditionalOptions";
	
	protected final static int defaultThreads = 4;
	
	
	private final SettingsModelString m_fastqc = 
			new SettingsModelString(CFGKEY_FASTQC,"");
	
	private final SettingsModelString m_outfolder = 
			new SettingsModelString(CFGKEY_OUTFOLDER,"");
	
	private final SettingsModelIntegerBounded m_threads = 
			new SettingsModelIntegerBounded(CFGKEY_THREADS, defaultThreads, 1, Integer.MAX_VALUE);
	
	private final SettingsModelString m_additional_options = 
			new SettingsModelString(CFGKEY_ADDITIONAL_OPTIONS,"");
	
	//The Output Col Names
	public static final String OUT_READFILE1 = "Path2ReadFile1";
	public static final String OUT_READFILE2 = "Path2ReadFile2";
	public static final String OUT_FILTERFILE = "Path2filterfile";
	
	//ReadType: paired-end or single-end
	private static String readType = "";
	
	//IBISKNIMENodesPlugin.FASTQC

	/**
     * Constructor for the node model.
     */
    protected FastQC_v2NodeModel() {
    
        super(1, 1);
        
        addSetting(m_fastqc);
    	addSetting(m_outfolder);
    	addSetting(m_threads);
    	addSetting(m_additional_options);
    	
    	addPrefPageSetting(m_fastqc, IBISKNIMENodesPlugin.FASTQC);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	String inFile1 = inData[0].iterator().next().getCell(0).toString();
    	String inFile2 = "";
    	
    	String outFile = inFile1;
    	String outFile2 = "";
    	
    	// Beautify outfolderPath string
    	String outFolderPath = m_outfolder.getStringValue();
    	outFolderPath = outFolderPath.trim();
    	if(!outFolderPath.equals("") && !outFolderPath.endsWith(System.getProperty("file.separator"))){
    		outFolderPath += System.getProperty("file.separator");
    	}
    	
    	int numThreads = m_threads.getIntValue();
    	
    	String additionalOpts = m_additional_options.getStringValue();
    	additionalOpts.trim();
    	
    	ArrayList<String> cmd = new ArrayList<String>();
    	cmd.add(IO.processFilePath(m_fastqc.getStringValue()));
    	
    	// Create outfile in outfolder if outfolder specified, otherwise at same location as infile
    	if(!outFolderPath.equals("")){
    		outFile = outFolderPath+inFile1.substring(inFile1.lastIndexOf(System.getProperty("file.separator"))+1, inFile1.length());
    		cmd.add("-o="+IO.processFilePath(outFolderPath));
    	} 
    	
    	
    	File lockFile = new File(outFile.substring(0,outFile.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	if(numThreads != 1){
    		cmd.add("-t="+numThreads);
    	}
    	
    	if(!additionalOpts.equals("")){
			setWarningMessage("NOTE! All additional FastQC options MUST use equals operator to bind arguments to options. i.e. --k=5");
			
			String[] addOpts = additionalOpts.split(" ");
    		for(String opt : addOpts){
    			cmd.add(opt);
    		}
		}
    	
    	cmd.add(inFile1);
    	
    	// Parse options into command string array - ensures spaces in path will be correctly interpreted
    	String[] command = new String[cmd.size()];
    	for(int indx=0; indx < cmd.size(); indx++){
    		command[indx] = cmd.get(indx);
    	}
    	
    	LOGGER.info("CMD: "+command.toString());
    	
    	LOGGER.info("Running FastQC...");
		LOGGER.info("Log files can be found in "+outFile+".stdOut and "+outFile+".stdErr");
		super.executeCommand(command, outFile, exec, lockFile, outFile+".stdOut", outFile+".stdErr");

		if(readType.equals("paired-end")){
			inFile2 = inData[0].iterator().next().getCell(1).toString();
			
			outFile2 = inFile2;
			
			// Create outfile in outfolder if outfolder specified, otherwise at same location as infile
	    	if(!outFolderPath.equals("")){
	    		outFile2 = outFolderPath+inFile2.substring(inFile2.lastIndexOf(System.getProperty("file.separator")), inFile2.length());
	    	} 
			
			lockFile = new File(outFile2.substring(0,outFile2.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);
			command[command.length-1] = inFile2;
			
			LOGGER.info("CMD: "+command.toString());
	    	
	    	LOGGER.info("Running FastQC...");
			LOGGER.info("Log files can be found in "+outFile2+".stdOut and "+outFile2+".stdErr");
			super.executeCommand(command, outFile2, exec, lockFile, outFile2+".stdOut", outFile2+".stdErr");
		}
		
		BufferedDataContainer cont;
	    FileCell[] c;
		
		cont = exec.createDataContainer(createSpecs());
	    
	    
	    if(readType.equals("single-end")){
	    	c = new FileCell[]{
		    		FileCellFactory.create(inFile1)};
	    } else {
	    	c = new FileCell[]{
	    			FileCellFactory.create(inFile1),
	       			FileCellFactory.create(inFile2)};
	    }
	    
	    cont.addRowToTable(new DefaultRow("Row0", c));
	  
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
		return new BufferedDataTable[]{outTable};
    }

	/**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	CompatibilityChecker CC = new CompatibilityChecker();
    	readType = CC.getReadType(inSpecs, 0);
    	if(CC.getWarningStatus()){
    		setWarningMessage(CC.getWarningMessages());
    	}
    	
    	super.updatePrefs();
        return new DataTableSpec[]{createSpecs()};
    }
    
    
    /**
     * Create Tablespecs
     * @return
     */
    private DataTableSpec createSpecs(){
    	DataTableSpec out;
    	
    	if(readType.equals("single-end")){
    		out = new DataTableSpec(
            		new DataColumnSpec[]{
            				new DataColumnSpecCreator(OUT_READFILE1, FastQCell.TYPE).createSpec()});
    	} else {
    		out = new DataTableSpec(
            		new DataColumnSpec[]{
            				new DataColumnSpecCreator(OUT_READFILE1, FastQCell.TYPE).createSpec(),
            				new DataColumnSpecCreator(OUT_READFILE2, FastQCell.TYPE).createSpec()});
    	}
    	
    	return out;
    }
 
}

