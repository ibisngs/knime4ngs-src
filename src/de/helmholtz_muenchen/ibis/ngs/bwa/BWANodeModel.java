package de.helmholtz_muenchen.ibis.ngs.bwa;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.workflow.NodeProgressEvent;
import org.knime.core.node.workflow.NodeProgressListener;

import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;



/**
 * This is the model implementation of BWA.
 * 
 */
public class BWANodeModel extends NodeModel implements NodeProgressListener {
    
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_BWAFILE = "bwafile";
	public static final String CFGKEY_CHECKCOLORSPACED = "checkColorSpaced";
	public static final String CFGKEY_BWTINDEX = "bwtIndex";
	public static final String CFGKEY_READTYPE = "readType";
	public static final String CFGKEY_CHECKINDEX = "checkIndexRefSeq";
	public static final String CFGKEY_READGROUP = "readgroup";
	public static final String CFGKEY_READGROUPBOOLEAN = "readgroupboolean";
	public static final String CFGKEY_ALNALGO = "alnalgo";
	
	private final SettingsModelString m_refseqfile = new SettingsModelString(BWANodeModel.CFGKEY_REFSEQFILE,"");
	private final SettingsModelString m_bwafile = new SettingsModelString(BWANodeModel.CFGKEY_BWAFILE,"");
	private final SettingsModelBoolean m_checkIndexRefSeq = new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKINDEX, true);
	private final SettingsModelBoolean m_checkColorSpaced = new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKCOLORSPACED	, false);
	private final SettingsModelString m_bwtIndex = new SettingsModelString(BWANodeModel.CFGKEY_BWTINDEX, "BWT-SW");
	private final SettingsModelString m_alnalgo = new SettingsModelString(BWANodeModel.CFGKEY_ALNALGO, "BWA-backtrack");
	private final SettingsModelString m_readType = new SettingsModelString(BWANodeModel.CFGKEY_READTYPE, "auto-detect");
	private final SettingsModelString m_readGroup = new SettingsModelString(CFGKEY_READGROUP, "");
	private final SettingsModelBoolean m_readGroupBoolean = new SettingsModelBoolean(CFGKEY_READGROUPBOOLEAN, false);
	
	
//TODO: implement alnignment algos....
    /**
     * Constructor for the node model.
     */
    protected BWANodeModel() {
    	
    	super(1, 1);
    	
    	m_readType.setEnabled(false);
    	
    }
//    void org.knime.core.node.NodeProgressMonitor.checkCanceled() throws CanceledExecutionException
    
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	    
    	/**Initialize logfile**/
    	String fle = inData[0].iterator().next().getCell(0).toString();
    	String logfile = fle.substring(0,fle.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("BWA"));
    	/**end initializing logfile**/
    	
    	String path2bwa = m_bwafile.getStringValue();
    	String path2refFile = m_refseqfile.getStringValue();
    	String path2readFile = inData[0].iterator().next().getCell(0).toString();
    	String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	String alnalgo = m_alnalgo.getStringValue();
    	String basePath = path2readFile.substring(0,path2readFile.lastIndexOf('/')+1);
    	String outBaseName1 = path2readFile.substring(path2readFile.lastIndexOf("/")+1,path2readFile.lastIndexOf("."));
    	String outBaseName = outBaseName1;
    	String outBaseName2 = outBaseName1;
    	if(path2readFile2.length() > 1 && !path2readFile2.equals("na")) {
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2readFile.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}

    	String out2Name = basePath+outBaseName+"_aln.sam";
    	String out1Name = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";
    	String colorSpaced = "-c ";


    	String samse = " samse ";
    	String sampe = " sampe ";
    	String bwasw = " bwasw ";
    	String bwamem = " bwamem ";
    	Boolean isBam = false;
    	
    	//Add Read Group if necessary:
    	if(m_readGroupBoolean.getBooleanValue()){
    		samse+="-r "+m_readGroup.getStringValue();
    		sampe+="-r "+m_readGroup.getStringValue();
    	}
    	
    	
    	if(path2readFile.substring(path2readFile.length()-3, path2readFile.length()) == "bam") {
    		isBam = true;
    	}
    	
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType)) {
    		readType = readTypePrevious;
    	}
    	logBuffer.append("Read Type: " + readType + "\n");
    	   	
    	/**
    	 * Create Index for Read Files
    	 */

    	bwa_index(exec,logBuffer, colorSpaced, path2bwa, path2refFile, path2readFile2);

    	
    	
    	/**
    	 * Run bwa aln
    	 */
     //    bwa aln sequence.fasta s_1_1_sequence.txt > aln_sa.sai
    	logBuffer.append("Find the SA coordinates of the input reads.\n");
    	bwa_aln(exec,readType, basePath, outBaseName, outBaseName1, outBaseName2, path2refFile, path2bwa, path2readFile, logBuffer, path2readFile2, isBam);
    	System.out.println("Finished BWA aln...");
    	
    	/**
    	 * Align reads
    	 */
    	
    	if(isBam) {
			path2readFile2 = path2readFile;
		}
// ### BWA-Backtrack ###
    	if(alnalgo.equals("BWA-backtrack")) {
    			bwa_backtrack(exec,readType,logBuffer,path2bwa,path2refFile,path2readFile,out1Name,out2Name,out11Name,out12Name,path2readFile2,samse,sampe);
// ### BWA-SW ###
    	} else if(alnalgo.equals("BWA-SW")) {
    			bwa_bwasw(exec,logBuffer, path2bwa, bwasw, path2refFile, path2readFile, path2readFile2, out2Name);
// ### BWA-MEM ###
    	} else if(alnalgo.equals("BWA-MEM")) {
    			bwa_mem(exec,logBuffer, path2bwa, bwamem, path2refFile, path2readFile, path2readFile2, out2Name);
    	}

    	DataColumnSpecCreator col = new DataColumnSpecCreator("Path2SAMFile", StringCell.TYPE);
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2RefFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col.createSpec(),col1.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(out2Name);
    	StringCell cl2 = new StringCell(path2refFile);
    	DataCell[] c = new DataCell[]{cl1,cl2};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();

    	pushFlowVariableString("BAMSAMINFILE",out2Name);
    	
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
     * Runs bwa index
     * @param logBuffer
     * @param colorSpaced
     * @param path2bwa
     * @param path2refFile
     * @param path2readFile
     * @throws Exception 
     */
    private void bwa_index(ExecutionContext exec,StringBuffer logBuffer, String colorSpaced, String path2bwa, String path2refFile, String path2readFile) throws Exception{
    	 // bwa index -a bwtsw sequence.fasta
    	
    	//Constant values
    	String indexBWTSW = " index -a bwtsw ";
    	String indexIS = " index -a is ";
    	String index;
    	String com = "";
    	
    	if(m_checkIndexRefSeq.getBooleanValue()) {
	    	logBuffer.append("Index reference sequence.\n");
	    	if(!m_checkColorSpaced.getBooleanValue()) {
	    		colorSpaced = "";
	    	}
	    	if(m_bwtIndex.getStringValue() == "BWT-SW") {
	    		index = indexBWTSW;
	    	} else {
	    		index = indexIS;
	    	}
	    	com = path2bwa + index + colorSpaced + path2refFile;
	    	// begin QueueSub #################################################
			if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
				String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWAindex";
				String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
				String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
				new QSub(com, name, memory, logfle, true);
	 			logBuffer.append("QSub: " + com + "\n");
				logBuffer.append("See external logfile: " + logfle + "\n");
			// end QueueSub ###################################################
			} else {
	    	//Execute
//	    		Executor.executeCommand(com, "BWA Indexing", exec);
			}
    	} else {
    		logBuffer.append("Indexing reference sequence SKIPPED.\n");
    	}
    }
    
    private void bwa_aln(ExecutionContext exec, String readType, String basePath, String outBaseName, String outBaseName1, String outBaseName2, String path2refFile, String path2bwa, String path2readFile, StringBuffer logBuffer, String path2readFile2, boolean isBam) throws Exception{
    	
    	String outName = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";
    	String com = ""; 
    	String outfile = outName;
    	
    	//System.out.println(path2readFile);
    	//System.out.println(path2readFile2);
    	
    	//If Inputfile is in bam format
    	String bam = "-b0 ";
    	String bam2 = "-b2 ";
    	if(!isBam) {
			bam = "";
			bam2 = "";
		}
    	if(readType.equals("paired-end")){
    		outfile = out11Name;
    		if(isBam){
    			bam = "-b1";
    		}
    	}
    	//Perform aln for forward reads OR single end reads
    	com = path2bwa + " aln "+bam+" "+ path2refFile +" "+ path2readFile +" -f "+ outfile;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWAaln";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	//Execute
//			Executor.executeCommand(com, "BWA Aln forward", exec);
			
		}
		//If paired end, repeat previous step
    	if(readType.equals("paired-end")) {
    		com = path2bwa + " aln "+bam2+" "+ path2refFile +" "+ path2readFile2 +" -f "+ out12Name;
    		// begin QueueSub #################################################
			if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
				String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWAaln";
				String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
				String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
				new QSub(com, name, memory, logfle, true);
	 			logBuffer.append("QSub: " + com + "\n");
				logBuffer.append("See external logfile: " + logfle + "\n");
			// end QueueSub ###################################################
			} else {
		    	//Execute
//	    		Executor.executeCommand(com, "BWA Aln Reverse", exec);
			}
    	}
    }
    
    
    private void bwa_backtrack(ExecutionContext exec,String readType, StringBuffer logBuffer, String path2bwa, String path2refFile, String path2readFile, String out1Name, String out2Name, String out11Name, String out12Name, String path2readFile2, String samse, String sampe) throws Exception{
		String na = "";
		String com = "";

    	if(readType.equals("single-end")) {
    // bwa samse sequence.fasta aln_sa.sai s_1_1_sequence.txt > aln.sam
    		logBuffer.append("Generate alignments in the SAM format given single-end reads.\n");
    		com = path2bwa + samse +" -f "+ out2Name +" "+ path2refFile +" "+ out1Name +" "+ path2readFile;
    		na = "samse";
    	} else {
    // bwa sampe sequence.fasta aln_sa_1.sai aln_sa_2.sai s_1_1_sequence.fq s_1_2_sequence.fq > aln.sam
    		logBuffer.append("Generate alignments in the SAM format given paired-end reads.\n");
    		com = path2bwa + sampe +" -f "+ out2Name +" "+ path2refFile +" "+ out11Name +" "+ out12Name +" "+ path2readFile +" "+ path2readFile2;
    		na = "sampe";
    	}
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWA" + na;
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	//Execute
//    		Executor.executeCommand(com, "BWA Backtrack", exec);
		}
    }
    
    private void bwa_bwasw(ExecutionContext exec, StringBuffer logBuffer, String path2bwa, String bwasw, String path2refFile, String path2readFile, String path2readFile2, String out2Name) throws Exception{
    	String com = "";
    	
		logBuffer.append("Generate alignments in the SAM format.\n");
		com = path2bwa + bwasw + " -f " + out2Name +" "+ path2refFile + " " + path2refFile + " " + path2readFile + " " + path2readFile2;
		// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWA-SW";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	//Execute
//    		Executor.executeCommand(com, "BWA SW", exec);
		}
    }
    
    private void bwa_mem(ExecutionContext exec, StringBuffer logBuffer, String path2bwa, String bwamem, String path2refFile, String path2readFile, String path2readFile2, String out2Name) throws Exception{
    	String com = "";
		logBuffer.append("Generate alignments in the SAM format.\n");
		com = path2bwa + bwamem + path2refFile + " " + path2refFile + " " + path2readFile + " " + path2readFile2;
		// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_BWA-MEM";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {	
	    	//Execute
//    		Executor.executeCommand(com, "BWA MEM", exec,out2Name);
		}
    }
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	// Warning if there is a problem with readType
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType) && !readType.equals("auto-detect")) {
    		setWarningMessage("The previous node indicates that you have " + readTypePrevious + " reads, but you have chosen " + readType + ". BWA will use " + readTypePrevious + " mapping.");
    	}
    	
        //Version control
        if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==1){
        	setWarningMessage("WARNING: You are using a newer BWA version than "+FileValidator.BWA_VERSION +"! This may cause problems");
        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==2){
        	throw new InvalidSettingsException("You are using a outdated version of BWA! Please update your version");
        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==-1){
        	System.out.println("Something wrong here");
        }
    	
    	if(m_refseqfile.getStringValue().length() > 1) {
	    	if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	            throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
	    	}
    	}
    	
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
    	m_bwafile.saveSettingsTo(settings);
    	m_refseqfile.saveSettingsTo(settings);
    	m_bwtIndex.saveSettingsTo(settings);
    	m_readType.saveSettingsTo(settings);
    	m_checkColorSpaced.saveSettingsTo(settings);
    	m_checkIndexRefSeq.saveSettingsTo(settings);
    	m_readGroup.saveSettingsTo(settings);
    	m_readGroupBoolean.saveSettingsTo(settings);
    	m_alnalgo.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_bwafile.loadSettingsFrom(settings);
    	m_refseqfile.loadSettingsFrom(settings);
    	m_bwtIndex.loadSettingsFrom(settings);
    	m_readType.loadSettingsFrom(settings);
    	m_checkColorSpaced.loadSettingsFrom(settings);
    	m_checkIndexRefSeq.loadSettingsFrom(settings);
    	m_readGroup.loadSettingsFrom(settings);
    	m_readGroupBoolean.loadSettingsFrom(settings);
    	m_alnalgo.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_bwafile.validateSettings(settings);
    	m_refseqfile.validateSettings(settings);
    	m_bwtIndex.validateSettings(settings);
    	m_readType.validateSettings(settings);
    	m_checkColorSpaced.validateSettings(settings);
    	m_checkIndexRefSeq.validateSettings(settings);
    	m_readGroup.validateSettings(settings);
    	m_readGroupBoolean.validateSettings(settings);
    	m_alnalgo.validateSettings(settings);
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



	@Override
	public void progressChanged(NodeProgressEvent pe) {
		
		
	}

}

