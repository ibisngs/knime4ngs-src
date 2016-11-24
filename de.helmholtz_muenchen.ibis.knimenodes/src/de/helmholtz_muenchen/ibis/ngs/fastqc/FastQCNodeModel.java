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
package de.helmholtz_muenchen.ibis.ngs.fastqc;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

/**
 * This is the model implementation of FastQC.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class FastQCNodeModel extends HTExecutorNodeModel {
    
	
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastQCNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2ReadFile1";
	public static final String OUT_COL2 = "Path2ReadFile2";
	public static final String OUT_COL3 = "Path2filterfile";
	
	//ReadType: paired-end or single-end
	private static String readType = "";
	
    /**
     * Constructor for the node model.
     */
    protected FastQCNodeModel() {
    
        super(1, 1);
        
    }

    /**
     * {@inheritDoc}
     */
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	   	
    	/** Get the input columns **/
    	String readsFile1 = inData[0].iterator().next().getCell(0).toString();
    	String outfile1 = this.getSettingsFileName(readsFile1);
    	String outfileMerged = outfile1;  	
    	
    	File lockFile = new File(readsFile1.substring(0,readsFile1.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);

    	/**Prepare Command**/
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("java");
    	String jarCall 			= "-jar "+IO.getScriptPath()+"libs/FastQC.jar ";
    	String path2mergeScript = "sh "+IO.getScriptPath() + "scripts/bash/mergeFsettings.sh";
    	command.add(jarCall + readsFile1);
    	
    	/**Execute for first file**/
    	String[] com = command.toArray(new String[command.size()]);
    	StringBuffer sysErr = new StringBuffer(50);
    	
    	super.executeCommand(new String[]{StringUtils.join(com, " ")}, exec, null, lockFile, null, null, null, sysErr, null);
    	
    	//Show FastQC Output
    	LOGGER.info(sysErr);
    	
    	/**If Paired-End data**/
    	String readsFile2 = "";
    	if(readType.equals("paired-end")) {
    		readsFile2 = inData[0].iterator().next().getCell(1).toString();
    		if (!readsFile2.equals("") && !readsFile2.equals(readsFile1)) {
    		
    			String outfile2 = this.getSettingsFileName(readsFile2); // override this path
    			String readFile1Path = IO.removeZipExtension(readsFile1);
    			String readFile2Name = IO.removeZipExtension(new File(readsFile2).getName());
    			
    			outfileMerged = readFile1Path.substring(0,readFile1Path.lastIndexOf(".")) + "_" + readFile2Name.substring(0,readFile2Name.lastIndexOf(".")) + "_fastqc.filterSettings";
	    		
	    		//Set new lock file for reverse read
	    		lockFile = new File(readsFile2.substring(0,readsFile2.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);
	    		
	    		//Replace readsFile1 with readsFile2 and execute again
	    		com[com.length-1] = jarCall + readsFile2;
	    		//Clear StringBuffer
	    		sysErr.setLength(0);
	    		sysErr.append("\n");
		    	super.executeCommand(new String[]{StringUtils.join(com, " ")}, exec, null, lockFile, null, null, null, sysErr, null);
//	    		//Show FastQC Output
	        	LOGGER.info(sysErr);
	        	//Clear StringBuffer
	    		sysErr.setLength(0);
	    		sysErr.append("\n");
	    		
	    		/** merge the two filter settings files */
	        	ArrayList<String> commandMerge = new ArrayList<String>();
//	        	commandMerge.add("sh");
	        	commandMerge.add(path2mergeScript);
	        	commandMerge.add(outfile1);
	        	commandMerge.add(outfile2);
	        	commandMerge.add(outfileMerged);
	        	        
	        	//No HTE, Merge is performed each time => strange null pointer exception
//	        	Executor.executeCommand(new String[]{StringUtils.join(commandMerge, " ")},exec,LOGGER,sysErr,sysErr);
	        	
	        	//Set new lock file for merging
	        	lockFile = new File(outfileMerged+SuccessfulRunChecker.LOCK_ENDING);
	        	super.executeCommand(new String[]{StringUtils.join(commandMerge, " ")},exec,
	        			null, lockFile, null, null , null, sysErr, null);

	        	
//	        	//Show FastQC Output
	        	LOGGER.info(sysErr);
			}
    	}

	    BufferedDataContainer cont;
	    FileCell[] c;
	    	
	    if(readType.equals("single-end")){
	        cont = exec.createDataContainer(createSpecs());
	        c = new FileCell[]{
	        		FileCellFactory.create(readsFile1),
	       			FileCellFactory.create(outfileMerged)};
	    		
	    }else{
	       	cont = exec.createDataContainer(createSpecs());
	       	c = new FileCell[]{
	       			FileCellFactory.create(readsFile1),
	       			FileCellFactory.create(readsFile2),
	       			FileCellFactory.create(outfileMerged)};
	   	}
	  
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	/**Push FlowVars**/
//    	String secFile = "";
//    	if(readsFile1.substring(readsFile1.length()-3,readsFile1.length()).equals("bam")) {
//    		secFile = "true";
//    	}else {
//    		secFile = "false";
//    	}
//    	pushFlowVariableString("isBAM", secFile);
    	
    	/**Delete FastQC zip files**/
    	deleteZipFiles(readsFile1, readsFile2);

    	
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
    	
    	return new DataTableSpec[]{createSpecs()};
    }

    /**
     * Deletes the FASTQC zip files
     * @param readsFile1
     * @param readsFile2
     */
    private void deleteZipFiles(String readsFile1, String readsFile2){
    	readsFile1 = IO.removeZipExtension(readsFile1);
    	File zipFile = new File(readsFile1 + "_fastqc.zip");
    	if(zipFile.exists()) {
    		zipFile.delete();
    	}
    	File zipFile1 = new File(readsFile1.substring(0,readsFile1.lastIndexOf(".")) + "_fastqc.zip");
    	if(zipFile1.exists()) {
    		zipFile1.delete();
    	}
    	if(!readsFile2.equals("")) {
    		readsFile2 = IO.removeZipExtension(readsFile2);
    		File zipFile2 = new File(readsFile2 + "_fastqc.zip");
        	if(zipFile2.exists()) {
        		zipFile2.delete();
        	}
        	File zipFile3 = new File(readsFile2.substring(0,readsFile1.lastIndexOf(".")) + "_fastqc.zip");
        	if(zipFile3.exists()) {
        		zipFile3.delete();
        	}
    	}
    }
    
    private String getSettingsFileName(String infile) {
    	infile = IO.removeZipExtension(infile);
    	infile = infile.replaceAll("\\.fq$", ".fq_fastqc.filterSettings");
    	infile = infile.replaceAll("\\.fastq$", "_fastqc.filterSettings");
    	return infile;
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
        					new DataColumnSpecCreator(OUT_COL1, FastQCell.TYPE).createSpec(),
        					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()});
    	}else{
    		out = new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL1, FastQCell.TYPE).createSpec(),
        					new DataColumnSpecCreator(OUT_COL2, FastQCell.TYPE).createSpec(),
        					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()});
    	}
    	return out;
    }
}

