package de.helmholtz_muenchen.ibis.ngs.fastqc;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;




/**
 * This is the model implementation of FastQC.
 * 
 *
 * @author Max
 */
public class FastQCNodeModel extends NodeModel {
    
	
    /**
     * Constructor for the node model.
     */
    protected FastQCNodeModel() {
    
        super(1, 1);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String readsFile1 = inData[0].iterator().next().getCell(0).toString();
    	String readsFile2 = inData[0].iterator().next().getCell(1).toString();
    	String readType = getAvailableInputFlowVariables().get("readType").getStringValue();
    	
    	/**Initialize logfile**/
    	String logfile = readsFile1.substring(0,readsFile1.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("FastQC"));
    	/**end initializing logfile**/

    	String path = FastQCNodeModel.class.getProtectionDomain().getCodeSource().getLocation().getPath();
    	
    	/**
    	 * AB HIER NEU !
    	 */
    	String sub =path.substring(path.lastIndexOf("/")+1, path.length());
    	String com="";
    	if(sub.equals("")){
    		com="cd "+path+"libs/FastQC_rebuild/ ; sh run_fastq2.sh "+readsFile2;
    	}else{//From Jar
    		String tmpfolder = path.substring(0, path.lastIndexOf("/")+1);
    		//com = "cd "+tmpfolder+" ; unzip -n "+sub+" ; sh FastQC_rebuild/run_fastq2.sh "+readsFile2;  #Jan#
    		com = "cd "+tmpfolder+"libs/FastQC_rebuild/ ; sh run_fastq2.sh "+readsFile2;
    	}	
    	
     //ALT !---->   String com = "cd " + path + "libs/FastQC_rebuild/ ; sh run_fastq2.sh " + readsFile2;
	/**
	 * Bis HIER !
	 */
    	
    	if(readType.equals("paired-end") && !readsFile2.equals("") && !readsFile2.equals(readsFile1)) {
			if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
				String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_FastQC-2";
				String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
				String logfle = readsFile1.substring(0,readsFile1.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
				new QSub(com, name, memory, logfle, false);
	 			logBuffer.append("QSub: " + com + "\n");
				logBuffer.append("See external logfile: " + logfle + "\n");
			// end QueueSub ###################################################
			} else {
		     	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
		    	Process p1 = b.start();
		    	p1.waitFor();
		        logBuffer.append(ShowOutput.getLogEntry(p1,com));
			}
		}
		/**
		 * UND AB HIER
		 */
    	if(sub.equals("")){
    		com="cd "+path+"libs/FastQC_rebuild/ ; sh run_fastq2.sh "+readsFile1;
    	}else{//From Jar
    		String tmpfolder = path.substring(0, path.lastIndexOf("/")+1);
    		//com = "cd "+tmpfolder+" ; unzip -n "+sub+" ; sh FastQC_rebuild/run_fastq2.sh "+readsFile1;	#Jan#
    		com = "cd "+tmpfolder+"libs/FastQC_rebuild/ ; sh run_fastq2.sh "+readsFile1;
    	}	
    	
		//OLD-->   com = "cd " + path + "libs/FastQC_rebuild/ ; sh run_fastq2.sh " + readsFile1;
    	
    	/**
    	 * BIS HIER !
    	 */
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_FastQC-1";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = readsFile1.substring(0,readsFile1.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	     	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p1 = b.start();
	    	p1.waitFor();
	        logBuffer.append(ShowOutput.getLogEntry(p1,com));
		}
    	
        //Create Output
        String outfile = readsFile1 + "_fastqc.filterSettings";
        String outfile2 = readsFile1.substring(0,readsFile1.lastIndexOf(".")) + "_fastqc.filterSettings";
        if(new File(outfile2).exists()) {
        	outfile = outfile2;
        }
        
        DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2ReadFile1", StringCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2ReadFile2", StringCell.TYPE);
        DataColumnSpecCreator col3 = new DataColumnSpecCreator("Path2filterfile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec(),col3.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(readsFile1);
    	StringCell cl2 = new StringCell(readsFile2);
    	StringCell cl3 = new StringCell(outfile);
    	DataCell[] c = new DataCell[]{cl1,cl2,cl3};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
    	String secFile = "";
    	//System.out.println(readsFile1.substring(readsFile1.length()-3,readsFile1.length()));
    	if(readsFile1.substring(readsFile1.length()-3,readsFile1.length()).equals("bam")) {
    		secFile = "true";
    	} else {
    		secFile = "false";
    	}
    	pushFlowVariableString("isBAM", secFile);
    	
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
    	
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
        return new BufferedDataTable[]{out};
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
    	
    	// Check input ports
    	String[] cn=inSpecs[0].getColumnNames();
    	if(!cn[0].equals("") && !cn[0].equals("Path2ReadFile1")) {
    		throw new InvalidSettingsException("This node is incompatible with the previous node. The outport of the previous node has to fit to the inport of this node.");
    	}
    	
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
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

