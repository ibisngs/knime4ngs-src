package de.helmholtz_muenchen.ibis.ngs.fastqc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;




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
    	String readsFile2 = inData[0].iterator().next().getCell(1).toString();
    	  
    	
//    	String readType = getAvailableInputFlowVariables().get("readType").getStringValue();
    	
        /**Create Output Specs**/
        String outfile1 = readsFile1.substring(0,readsFile1.lastIndexOf(".")) + "_fastqc.filterSettings";
        
        /**
         * fq format handling. fastq is already handled in Fastqc.jar
         */
        if (readsFile1.substring(readsFile1.length()-2, readsFile1.length()).equals("fq"))
    	{
        	outfile1 = readsFile1.substring(0,readsFile1.lastIndexOf(".")) + ".fq_fastqc.filterSettings";
    	}
        
        String outFileSettings = outfile1;
    	String outfileMerged = outfile1;
    	String readFile2Name;
    	  	
    	
    	/**If Paired-End data**/
    	if(readType.equals("paired-end") && !readsFile2.equals("") && !readsFile2.equals(readsFile1)) {
    		readFile2Name = new File(readsFile2).getName();
        	outFileSettings = outfileMerged; // override this path
        	outfileMerged = readsFile1.substring(0,readsFile1.lastIndexOf(".")) + "_" + readFile2Name.substring(0,readFile2Name.lastIndexOf(".")) + "_fastqc.filterSettings";
    	}
    	
    	File lockFile = new File(readsFile1.substring(0,readsFile1.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);

	  	
	    	/**Prepare Command**/
	    	ArrayList<String> command = new ArrayList<String>();
	    	command.add("java");
	    	String jarCall 			= "-jar "+IO.getScriptPath()+"libs/FastQC.jar ";
	    	String path2mergeScript = IO.getScriptPath() + "scripts/bash/mergeFsettings.sh";
	    	command.add(jarCall + readsFile1);
	
	    	
	    	/**Execute for first file**/
	    	String[] com = command.toArray(new String[command.size()]);
	    	StringBuffer sysErr = new StringBuffer(50);
	    	
	    	super.executeCommand(new String[]{StringUtils.join(com, " ")}, exec, null, lockFile, null, null, null, sysErr, null);
	    	
	    	//Show FastQC Output
	    	LOGGER.info(sysErr);

	    	/**If Paired-End data**/
	    	if(readType.equals("paired-end") && !readsFile2.equals("") && !readsFile2.equals(readsFile1)) {
	    		
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
	        	String outfile2 = readsFile2.substring(0,readsFile2.lastIndexOf(".")) + "_fastqc.filterSettings";
	        	if (readsFile2.substring(readsFile2.length()-2, readsFile2.length()).equals("fq"))
	         	{
	             	outfile2 = readsFile2.substring(0,readsFile2.lastIndexOf(".")) + ".fq_fastqc.filterSettings";
	         	}
	        	// merge the two settings files
	        	ArrayList<String> commandMerge = new ArrayList<String>();
	        	commandMerge.add("sh");
	        	commandMerge.add(path2mergeScript);
	        	commandMerge.add(outfile1);
	        	commandMerge.add(outfile2);
	        	commandMerge.add(outfileMerged);
	        	        
	        	//No HTE, Merge is performed each time
	        	Executor.executeCommand(new String[]{StringUtils.join(commandMerge, " ")},exec,LOGGER,sysErr,sysErr);
	        	
//	        	//Show FastQC Output
	        	LOGGER.info(sysErr);
			}

	    	BufferedDataContainer cont;
	    	FileCell[] c;
	    	
	    	if(readType.equals("single-end")){
	        	cont = exec.createDataContainer(createSpecs());
	        	c = new FileCell[]{
	        			(FileCell) FileCellFactory.create(readsFile1),
	        			(FileCell) FileCellFactory.create(outFileSettings)};
	    		
	    	}else{
	        	cont = exec.createDataContainer(createSpecs());
	        	c = new FileCell[]{
	        			(FileCell) FileCellFactory.create(readsFile1),
	        			(FileCell) FileCellFactory.create(readsFile2),
	        			(FileCell) FileCellFactory.create(outFileSettings)};
	    	}
	  
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	/**Push FlowVars**/
    	String secFile = "";
    	if(readsFile1.substring(readsFile1.length()-3,readsFile1.length()).equals("bam")) {
    		secFile = "true";
    	}else {
    		secFile = "false";
    	}
    	pushFlowVariableString("isBAM", secFile);
    	
    	/**Delete FastQC zip files**/
    	deleteZipFiles(readsFile1, readsFile2);

    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	
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
    	File zipFile = new File(readsFile1 + "_fastqc.zip");
    	if(zipFile.exists()) {
    		zipFile.delete();
    	}
    	File zipFile1 = new File(readsFile1.substring(0,readsFile1.lastIndexOf(".")) + "_fastqc.zip");
    	if(zipFile1.exists()) {
    		zipFile1.delete();
    	}
    	if(!readsFile2.equals("")) {
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
    
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }
}

