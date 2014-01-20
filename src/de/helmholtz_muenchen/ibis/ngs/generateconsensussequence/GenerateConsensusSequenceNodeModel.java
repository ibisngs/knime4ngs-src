package de.helmholtz_muenchen.ibis.ngs.generateconsensussequence;

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

import de.helmholtz_muenchen.ibis.ngs.fastqc.FastQCNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;

/**
 * This is the model implementation of GenerateConsensusSequence.
 * 
 */
public class GenerateConsensusSequenceNodeModel extends NodeModel {
    
    /**
     * Constructor for the node model.
     */
    protected GenerateConsensusSequenceNodeModel() {
    
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String path2samtools = inData[0].iterator().next().getCell(0).toString();
    	String path2bamFile = inData[0].iterator().next().getCell(1).toString();
    	String path2seqFile = inData[0].iterator().next().getCell(2).toString();
    	String path2bcftools = path2samtools.substring(0,path2samtools.length()-8)+"bcftools/bcftools";
    	String path2vcfutils = path2samtools.substring(0,path2samtools.length()-8)+"bcftools/vcfutils.pl";
    	String path2consensusfilefq = path2bamFile.substring(0,path2bamFile.lastIndexOf("."))+"_consensusSequence.fq";
    	String path2consensusfilefa = path2bamFile.substring(0,path2bamFile.lastIndexOf("."))+"_consensusSequence.fa";
    	String basepath = path2bamFile.substring(0, path2bamFile.lastIndexOf("/"));
    	
    	/**Initialize logfile**/
    	String logfile = path2seqFile.substring(0,path2seqFile.lastIndexOf(".")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("GenerateConsensusSequence"));
    	/**end initializing logfile**/
    	
    // samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq
    	String com = path2samtools + " mpileup -uf " + path2seqFile + " " + path2bamFile + " | " + path2bcftools + " view -cg - | " + path2vcfutils + " vcf2fq > " + path2consensusfilefq;
    	System.out.println(com);
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_GCS-Mpileup";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamFile.substring(0,path2bamFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p = b.start();
	    	p.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p, com));
		}
    	
    // Converting fastq to fasta
    	String path = FastQCNodeModel.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        String com1 = "cd " + basepath + "; perl " + path + "libs/fastq2fasta.pl -a " + path2consensusfilefq;
        System.out.println(com1);
     // begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_GCS-fQ2fA";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamFile.substring(0,path2bamFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com1, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com1 + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	ProcessBuilder b1 = new ProcessBuilder("/bin/sh", "-c", com1);
	    	Process p1 = b1.start();
	    	p1.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p1, com1));
		}
    	
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2ConsensusFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(path2consensusfilefa);
    	DataCell[] c = new DataCell[]{cl1};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
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

