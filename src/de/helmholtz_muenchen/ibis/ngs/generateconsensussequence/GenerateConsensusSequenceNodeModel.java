package de.helmholtz_muenchen.ibis.ngs.generateconsensussequence;

import java.io.File;
import java.io.IOException;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.ngs.fastqc.FastQCNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of GenerateConsensusSequence.
 * 
 */
public class GenerateConsensusSequenceNodeModel extends NodeModel {
    
	//The Output Col Names
	public static final String OUT_COL1 = "Path2ConsensusFile";
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GenerateConsensusSequenceNodeModel.class);
	
	
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
//    	String basepath = path2bamFile.substring(0, path2bamFile.lastIndexOf("/"));
    	
    	/**Initialize logfile**/
    	String logfile = path2seqFile.substring(0,path2seqFile.lastIndexOf(".")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("GenerateConsensusSequence"));
    	/**end initializing logfile**/
    	
    // samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq
    	String com = "sh -c \"."+path2samtools + " mpileup -uf " + path2seqFile + " " + path2bamFile + " | " + path2bcftools + " view -cg - | " + path2vcfutils + " vcf2fq > " + path2consensusfilefq+"\"";

    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(com, " ")},exec,LOGGER);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
    // Converting fastq to fasta
    	String path = FastQCNodeModel.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        String com1 = "perl " + path + "libs/fastq2fasta.pl -a " + path2consensusfilefq;

    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(com1, " ")},exec,LOGGER);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2consensusfilefa)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
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

