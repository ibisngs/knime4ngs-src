package de.helmholtz_muenchen.ibis.ngs.segemehl;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;

/**
 * This is the model implementation of Segemehl.
 * 
 */
public class SegemehlNodeModel extends HTExecutorNodeModel {
	
	public static final String CFGKEY_SEGEMEHLFILE = "segemehlfile";
	public static final String CFGKEY_ACCURACY= "segemehAccuracy";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_READTYPE = "readtype";
	public static final String CFGKEY_THREADS = "threads";
	public static final String CFGKEY_CLIP5ADAPTER = "clip5adapter";
	public static final String CFGKEY_CLIP3ADAPTER = "clip3adapter";
	public static final String CFGKEY_ADAPTER5SEQ = "adapter5seq";
	public static final String CFGKEY_ADAPTER3SEQ = "adapter3seq";
	public static final String CFGKEY_AUTOADAPTER3SEQ = "autoadapter3seq";
	public static final String CFGKEY_CLIPPOLYA = "clippolya";
	public static final String CFGKEY_CLIPPINGACCURACY = "clippingaccuracy";
	public static final String CFGKEY_SOFTHARDCLIPPING = "softhardclipping";
	public static final String CFGKEY_CHECKINDEX = "checkIndexRefSeq";
	public static final String CFGKEY_CHECKSPLITREADMAPPING = "checkSplitReadMapping";
	public static final String CFGKEY_CHECKSBISULFITEMAPPING = "checkBisulfiteMapping";
	public static final String CFGKEY_BISULFITEMAPPINGTYPE = "bisulfiteMappingType";

	private final SettingsModelString m_segemehlfile = new SettingsModelString(SegemehlNodeModel.CFGKEY_SEGEMEHLFILE,"");
	private final SettingsModelString m_refseqfile = new SettingsModelString(SegemehlNodeModel.CFGKEY_REFSEQFILE,"");
	private final SettingsModelString m_readType = new SettingsModelString(SegemehlNodeModel.CFGKEY_READTYPE,"");
	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_THREADS, 4, 1, 250);
	private final SettingsModelBoolean m_clip5adapter = new SettingsModelBoolean(CFGKEY_CLIP5ADAPTER, false);
	private final SettingsModelBoolean m_clip3adapter = new SettingsModelBoolean(CFGKEY_CLIP3ADAPTER, false);
	private final SettingsModelBoolean m_autoadapter3seq = new SettingsModelBoolean(CFGKEY_AUTOADAPTER3SEQ, false);
	private final SettingsModelString m_adapter5seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER5SEQ,"");
	private final SettingsModelString m_adapter3seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER3SEQ,"");
	private final SettingsModelBoolean m_clippolya = new SettingsModelBoolean(CFGKEY_CLIPPOLYA, false);
	private final SettingsModelIntegerBounded m_clippingaccuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_CLIPPINGACCURACY, 70, 0, 100);
	private final SettingsModelString m_softhardclipping = new SettingsModelString(SegemehlNodeModel.CFGKEY_SOFTHARDCLIPPING,"");
	private final SettingsModelBoolean m_checkIndexRefSeq = new SettingsModelBoolean(CFGKEY_CHECKINDEX, true);
	private final SettingsModelBoolean m_checkSplitReadMapping = new SettingsModelBoolean(CFGKEY_CHECKSPLITREADMAPPING, false);
	private final SettingsModelBoolean m_checkBisulfiteMapping = new SettingsModelBoolean(CFGKEY_CHECKSBISULFITEMAPPING, false);
	private final SettingsModelString m_bisulfiteMappingType = new SettingsModelString(SegemehlNodeModel.CFGKEY_BISULFITEMAPPINGTYPE,"");
	private final SettingsModelIntegerBounded m_accuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_ACCURACY, 90, 0, 100);

	//The Output Col Names
	public static final String OUT_COL1 = "Path2SAMFile";
	public static final String OUT_COL2 = "Path2RefFile";
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(SegemehlNodeModel.class);
		
    /**
     * Constructor for the node model.
     */
    protected SegemehlNodeModel() {
    
        super(1, 1);
        
        m_autoadapter3seq.setEnabled(false);
        m_adapter3seq.setEnabled(false);
        m_adapter5seq.setEnabled(false);
        m_clippingaccuracy.setEnabled(false);
        m_softhardclipping.setEnabled(false);
        m_bisulfiteMappingType.setEnabled(false);
        m_checkBisulfiteMapping.setEnabled(false);
        m_readType.setStringValue("single-end");
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	String path2reads1 = inData[0].iterator().next().getCell(0).toString();
    	
    	String path2segemehl = m_segemehlfile.getStringValue();
    	String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	String path2refSeq = m_refseqfile.getStringValue();
    	String path2indexedRefSeq = path2refSeq.substring(0,path2refSeq.lastIndexOf(".")+1)+"idx";
    	String basePath = path2reads1.substring(0,path2reads1.lastIndexOf('/')+1);
    	String outBaseName1 = path2reads1.substring(path2reads1.lastIndexOf("/")+1,path2reads1.lastIndexOf("."));
    	String outBaseName = outBaseName1;
    	String outBaseName2 = outBaseName1;
    	if(path2readFile2.length() > 1) {
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2reads1.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}
    	String outName = basePath+outBaseName+"_map.sam";
    	String outNameUnmatchedReads = path2refSeq.substring(0,path2refSeq.lastIndexOf("/")+1)+"unmatchedReads.f";
    	int nrOfThreads = m_threads.getIntValue();
    	int accuracy = m_accuracy.getIntValue();
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType)) {
    		readType = readTypePrevious;
    	}
    	
    // Indexing reference sequence: segemehl -x chr1.idx -d chr1.fa
    	if(m_checkIndexRefSeq.getBooleanValue()) {
    		LOGGER.info("Indexing reference sequence.");
    		command.add(path2segemehl);
    		command.add("-x "+path2indexedRefSeq);
	    	command.add("-d "+path2refSeq);
	     	/**
	     	 * Execute
	     	 */
	     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
	     	command = new ArrayList<String>();	//Clear Array
	    	
    	} else {
    		LOGGER.info("Indexing reference sequence SKIPPED.");
    	}
    	
    	
    // Match/ align reads: segemehl.x -i chr1.idx -d chr1.fa -q myreads.fa > mymap.sam
    	LOGGER.info("Match/ align reads to indexed reference sequence.");
    	
    	command.add(path2segemehl);
    	command.add("-s");
    	
    	if(m_clip5adapter.getBooleanValue()) {
    		command.add("-P " + m_adapter5seq.getStringValue());
    	}
    	if(m_clip3adapter.getBooleanValue()) {
    		if(m_autoadapter3seq.getBooleanValue()) {
    			command.add("-Y");
    		} else {
    			command.add("-Q " + m_adapter3seq.getStringValue());
    		}
    	}
    	if(m_clippolya.getBooleanValue()) {
    		command.add("-T");
    	}
    	if(m_clip5adapter.getBooleanValue() || m_clip3adapter.getBooleanValue() || m_clippolya.getBooleanValue()) {
    		command.add("-R " + m_clippingaccuracy.getIntValue());
    		if(m_softhardclipping.getStringValue().equals("Hard")) {
    			command.add("-C");
    		}
    	}
    	if(readType.equals("paired-end")) {
    		command.add("-p " + path2readFile2);
    	}
    	if(m_checkSplitReadMapping.getBooleanValue()) {
    		command.add("-S");
    	}
    	if(m_checkBisulfiteMapping.getBooleanValue()) {
    		String bmType = m_bisulfiteMappingType.getStringValue();
    		if(bmType.equals("methylC-seq/Lister et al.")) {
    			command.add("-F 1");
    		} else if(bmType.equals("bs-seq/Cokus et al. protocol")) {
    			command.add("-F 2");
    		} else if(bmType.equals("PAR-CLIP with 4SU")) {
    			command.add("-F 3");
    		} else if(bmType.equals("PAR-CLIP with 6SG")) {
    			command.add("-F 4");
    		}
    	}

    	command.add("-A "+accuracy);
    	command.add("-t "+nrOfThreads);
    	command.add("-i "+path2indexedRefSeq);
    	command.add("-d "+path2refSeq);
    	command.add("-q "+path2reads1);
    	command.add("-o "+outName);
    	command.add("-u "+outNameUnmatchedReads);
    	
     	/**
     	 * Execute
     	 */
    	File lockFile = new File(outName + ".klock");
	
		// execute the command
//		Executor.executeCommand(new String[]{joinedCommand},exec,LOGGER);
		super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, lockFile);


     	command = new ArrayList<String>();	//Clear Array
     	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outName),
    			(FileCell) FileCellFactory.create(path2refSeq)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	pushFlowVariableString("BAMSAMINFILE",outName);
    	
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
    	
    	// Warning if there is a problem with readType
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType)) {
    		setWarningMessage("The previous node indicates that you have " + readTypePrevious + " reads, but you have chosen " + readType + ". Segemehl will use " + readTypePrevious + " mapping.");
    	}
    	
    	//Warning concerning files
    	String isBam = getAvailableInputFlowVariables().get("isBAM").getStringValue();
    	if(isBam.equals("true")) {
    		if(readTypePrevious.equals("single-end")) {
    			throw new InvalidSettingsException("Segemehl does not support BAM files. Please choose a FastQ or FastA file containing your reads.");
    		} else {
    			throw new InvalidSettingsException("Segemehl does not support BAM files. Please choose two FastQ or FastA files containing your reads.");
    		}
    	}
    	
		if(!m_checkIndexRefSeq.getBooleanValue()) {
	    	String path2refSeq = m_refseqfile.getStringValue();
			File file = new File(path2refSeq.substring(0,path2refSeq.lastIndexOf(".")+1)+"idx");    	        
			if(!file.exists()){
				throw new InvalidSettingsException("The reference sequence has not been indexed yet. Please modify the options of Segemehl to run Segemehl.");
			}
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
    	/** added for HTE **/
    	super.saveSettingsTo(settings);
    	
    	m_segemehlfile.saveSettingsTo(settings);
    	m_refseqfile.saveSettingsTo(settings);
    	m_adapter3seq.saveSettingsTo(settings);
    	m_adapter5seq.saveSettingsTo(settings);
    	m_autoadapter3seq.saveSettingsTo(settings);
    	m_clip3adapter.saveSettingsTo(settings);
    	m_clip5adapter.saveSettingsTo(settings);
    	m_clippolya.saveSettingsTo(settings);
    	m_readType.saveSettingsTo(settings);
    	m_softhardclipping.saveSettingsTo(settings);
    	m_threads.saveSettingsTo(settings);
    	m_clippingaccuracy.saveSettingsTo(settings);
    	m_checkIndexRefSeq.saveSettingsTo(settings);
    	m_checkSplitReadMapping.saveSettingsTo(settings);
    	m_checkBisulfiteMapping.saveSettingsTo(settings);
    	m_bisulfiteMappingType.saveSettingsTo(settings);
    	m_accuracy.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	/** added for HTE **/
    	super.loadValidatedSettingsFrom(settings);
    	
    	m_segemehlfile.loadSettingsFrom(settings);
    	m_refseqfile.loadSettingsFrom(settings);
    	m_adapter3seq.loadSettingsFrom(settings);
    	m_adapter5seq.loadSettingsFrom(settings);
    	m_autoadapter3seq.loadSettingsFrom(settings);
    	m_clip3adapter.loadSettingsFrom(settings);
    	m_clip5adapter.loadSettingsFrom(settings);
    	m_clippolya.loadSettingsFrom(settings);
    	m_readType.loadSettingsFrom(settings);
    	m_softhardclipping.loadSettingsFrom(settings);
    	m_threads.loadSettingsFrom(settings);
    	m_clippingaccuracy.loadSettingsFrom(settings);
    	m_checkIndexRefSeq.loadSettingsFrom(settings);
    	m_checkSplitReadMapping.loadSettingsFrom(settings);
    	m_checkBisulfiteMapping.loadSettingsFrom(settings);
    	m_bisulfiteMappingType.loadSettingsFrom(settings);
    	m_accuracy.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	/** added for HTE **/
    	super.validateSettings(settings);
    	
    	m_segemehlfile.validateSettings(settings);
    	m_refseqfile.validateSettings(settings);
    	m_adapter3seq.validateSettings(settings);
    	m_adapter5seq.validateSettings(settings);
    	m_autoadapter3seq.validateSettings(settings);
    	m_clip3adapter.validateSettings(settings);
    	m_clip5adapter.validateSettings(settings);
    	m_clippolya.validateSettings(settings);
    	m_readType.validateSettings(settings);
    	m_softhardclipping.validateSettings(settings);
    	m_threads.validateSettings(settings);
    	m_clippingaccuracy.validateSettings(settings);
    	m_checkIndexRefSeq.validateSettings(settings);
    	m_checkSplitReadMapping.validateSettings(settings);
    	m_checkBisulfiteMapping.validateSettings(settings);
    	m_bisulfiteMappingType.validateSettings(settings);
    	m_accuracy.validateSettings(settings);
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

