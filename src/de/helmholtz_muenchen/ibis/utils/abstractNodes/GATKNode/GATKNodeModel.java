package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;

public abstract class GATKNodeModel extends HTExecutorNodeModel{

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_GATK_MEM = "GATK_MEM";
	public static final String CFGKEY_PATH2BED = "bed";
	public static final String CFGKEY_BED_FILE_CHECKBOX = "BED_FILE_CHECKBOX";
	public static final String CFGKEY_OPT_FLAGS ="OPT_FLAGS";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_path2bed = new SettingsModelString(GATKNodeModel.CFGKEY_PATH2BED,"");
    private final SettingsModelBoolean m_bed_file_checkbox = new SettingsModelBoolean(GATKNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
    private final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(CFGKEY_OPT_FLAGS,"",false);

	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "OUTFILE";
	
	private boolean outtable = true;
    
    /**
     * Constructor for the node model.
     */
    protected GATKNodeModel(PortType[] INPORTS, PortType[] OUTPORTS) {
        super(INPORTS, OUTPORTS);
        if(OUTPORTS.length==0) {
        	outtable = false;
        }
        
        addSetting(m_GATK);
   	 	addSetting(m_REF_GENOME);
   	 	addSetting(m_GATK_MEM);
   	 	addSetting(m_path2bed);
   	 	addSetting(m_bed_file_checkbox);
   	 	addSetting(m_OPT_FLAGS);
        
   	 	m_path2bed.setEnabled(false);
    }
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-Xmx"+m_GATK_MEM.getIntValue()+"G");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T "+getCommandWalker());
    	command.add("-R "+m_REF_GENOME.getStringValue());    	
    	command.add(getCommandParameters(inData)); //must be called before getOutfile() !!
   
    	String OUTFILE = getOutfile(); 	
    	command.add("-o "+OUTFILE);
    	
    	File lockFile = new File(OUTFILE+SuccessfulRunChecker.LOCK_ENDING); 
    	
    	if(m_bed_file_checkbox.getBooleanValue()){
			if(m_path2bed.getStringValue().equals("") || Files.notExists(Paths.get(m_path2bed.getStringValue()))) {
				setWarningMessage("Specify valid bed file!");
			} else {
				command.add("-L "+m_path2bed.getStringValue());
			}
    	}
    	
    	command.add(" "+m_OPT_FLAGS.getStringValue());
    	
    	LOGGER.info(StringUtils.join(command, " "));
     	super.executeCommand(new String[]{StringUtils.join(command, " ")},exec, lockFile,OUTFILE+".stdOut",OUTFILE+".stdErr");
     	
     	
    	/**
    	 * OUTPUT
    	 */
     	if(!outtable) {
     		return null;
     	}
     	
     	//Table1
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1_TABLE1, getOutColType()).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable1 = cont.getTable();


        return new BufferedDataTable[]{outTable1};
 
    }

    
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	if(!checkInputCellType(inSpecs)) {
    		new InvalidSettingsException("This node seems to be incompatible with the precedent node!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_GATK.getStringValue())) {
    		new InvalidSettingsException("Set path to GenomeAnalysisTK.jar!");
    	}

    	//check reference genome
    	String reffile = m_REF_GENOME.getStringValue();
    	if(CompatibilityChecker.inputFileNotOk(reffile)) {
    		new InvalidSettingsException("Set path to reference genome!");
    	}
    	
    	if (!Files.exists(Paths.get(reffile + ".fai"))) {
			throw new InvalidSettingsException("Reference sequence index: " + reffile + ".fai does not exist!");
		}

		String refbase = PathProcessor.getBase(reffile);
		if (!Files.exists(Paths.get(refbase + ".dict"))) {
			throw new InvalidSettingsException("Reference sequence dictionary: " + refbase + ".dict does not exist!");
		}
    	
		if (!outtable) {
			return null;
		}
		
		extraConfig();

		DataTableSpec outSpecTable1 = new DataTableSpec(
				new DataColumnSpec[] { new DataColumnSpecCreator(
						OUT_COL1_TABLE1, getOutColType()).createSpec() });
		return new DataTableSpec[] { outSpecTable1 };
    }
    

    /****************************** ABSTRACT METHODS **********************************/
    /**
     * Provides the node specific filter settings
     * @return
     */
    protected abstract String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException;
    	
    /**
     * Provides the GATK Walker
     * @return
     */
    protected abstract String getCommandWalker();
    protected abstract String getOutfile();
    protected abstract boolean checkInputCellType(DataTableSpec[] inSpecs);
    protected abstract DataType getOutColType();
    protected abstract void extraConfig() throws InvalidSettingsException;
}