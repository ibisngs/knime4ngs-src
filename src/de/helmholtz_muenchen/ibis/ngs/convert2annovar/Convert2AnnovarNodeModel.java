package de.helmholtz_muenchen.ibis.ngs.convert2annovar;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of Convert2Annovar.
 * 
 *
 * @author Sebastian Kopetzky
 */
public class Convert2AnnovarNodeModel extends NodeModel {
	
	public static final String CFGKEY_INSTALLPATH = "installpath";
	public static final String CFGKEY_INFILE = "infile";
	public static final String CFGKEY_SNPQUAL = "snpqual";
	public static final String CFGKEY_SNPPVALUE = "snppvalue";
	public static final String CFGKEY_COVERAGE = "coverage";
	public static final String CFGKEY_IFMAXCOVERAGE = "ifmaxcoverage";
	public static final String CFGKEY_MAXCOVERAGE = "maxcoverage";
	public static final String CFGKEY_INCLUDEINFO = "includeinfo";
	public static final String CFGKEY_ALLELICFRAC = "allelicfrac";
	public static final String CFGKEY_SPECIES = "species";
	public static final String CFGKEY_ALLALLELE = "allallele";
	public static final String CFGKEY_WITHZYG = "withzyg";

	public static final int DEFAULT_SNPQUAL = 20;
	public static final int DEFAULT_SNPPVALUE = 1;
	public static final int DEFAULT_COVERAGE = 0;
	public static final int DEFAULT_MAXCOVERAGE = 0;

	private final SettingsModelString m_installpath = new SettingsModelString(CFGKEY_INSTALLPATH,"");
	private final SettingsModelString m_infile = new SettingsModelString(CFGKEY_INFILE,"");
	private final SettingsModelIntegerBounded m_snpqual = new SettingsModelIntegerBounded(CFGKEY_SNPQUAL,DEFAULT_SNPQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_snppvalue = new SettingsModelDoubleBounded(CFGKEY_SNPPVALUE,DEFAULT_SNPPVALUE,0,1);
	private final SettingsModelIntegerBounded m_coverage = new SettingsModelIntegerBounded(CFGKEY_COVERAGE,DEFAULT_COVERAGE,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_ifmaxcoverage = new SettingsModelBoolean(CFGKEY_IFMAXCOVERAGE, false);
	private final SettingsModelIntegerBounded m_maxcoverage = new SettingsModelIntegerBounded(CFGKEY_MAXCOVERAGE,DEFAULT_MAXCOVERAGE,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_includeinfo = new SettingsModelBoolean(CFGKEY_INCLUDEINFO, false);
	private final SettingsModelBoolean m_allelicfrac = new SettingsModelBoolean(CFGKEY_ALLELICFRAC, false);
	private final SettingsModelBoolean m_species = new SettingsModelBoolean(CFGKEY_SPECIES, true);
	private final SettingsModelBoolean m_allallele = new SettingsModelBoolean(CFGKEY_ALLALLELE, false);
	private final SettingsModelBoolean m_withzyg = new SettingsModelBoolean(CFGKEY_WITHZYG, false);

	private static final NodeLogger LOGGER = NodeLogger.getLogger(Convert2AnnovarNodeModel.class);
	//The Output Col Names
	public static final String OUT_COL1 = "Path2OutFile";
	public static final String OUT_COL2 = "InstallPath";
	
	
    /**
     * Constructor for the node model.
     */
    protected Convert2AnnovarNodeModel() {
    	
        super(0, 1);
        
        m_snpqual.setEnabled(false);
        m_coverage.setEnabled(false);
        m_allelicfrac.setEnabled(false);
        m_species.setEnabled(false);
        m_snppvalue.setEnabled(false);
        m_allallele.setEnabled(false);
        m_withzyg.setEnabled(false);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String installPath = m_installpath.getStringValue();
    	String path2inFile = m_infile.getStringValue();
    	
    	/**Initialize logfile**/
    	String logfile = path2inFile.substring(0,path2inFile.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("BAMLoader"));
    	/**end initializing logfile**/
    	
    	ArrayList<String> command = new ArrayList<String>();

    	String fileExtension = path2inFile.substring(path2inFile.lastIndexOf(".")+1,path2inFile.length());
    	String path2outfile = path2inFile.substring(0,path2inFile.lastIndexOf(".")) + ".query";
    	String outfile = "-outfile " + path2outfile;
    	
    	command.add(installPath + "/convert2annovar.pl");
    	command.add("--verbose");
    	
    	if(fileExtension.equals("vcf")) {
    		command.add("-format vcf4");
    	} else if(fileExtension.equals("pileup")) {
    		command.add("-format pileup");
    	} else if(fileExtension.equals("tsv")) {
    		command.add("-format cg");
    	} else if(fileExtension.equals("gff")) {
    		command.add("-format gff3-solid");
    	}
    //	if(m_snpqual.isEnabled()) {
    //		snpqval = "-snpqual " + m_snpqual.getIntValue() + " ";
    //	}
    	if(m_snppvalue.isEnabled()) {
    		command.add("-snppvalue " + m_snppvalue.getDoubleValue());
    	}
    	if(m_coverage.isEnabled()) {
    		command.add("-coverage " + m_coverage.getIntValue());
    	}
    	if(m_maxcoverage.isEnabled()) {
    		command.add("-maxcoverage " + m_maxcoverage.getIntValue());
    	}
    	if(m_includeinfo.getBooleanValue()) {
    		command.add("-includeinfo");
    	}
    	if(m_allelicfrac.isEnabled()) {
    		if(m_allelicfrac.getBooleanValue()) {
    			command.add("-allelicfrac");
    		}
    	}
    	if(m_species.isEnabled()) {
    		if (m_species.getBooleanValue()) {
    			command.add("-species human");
    		}
    	}
    	if(m_allallele.isEnabled()) {
    		if(m_allallele.getBooleanValue()) {
    			command.add("-allallele");
    		}
    	}
    	if(m_withzyg.isEnabled()) {
    		if(m_withzyg.getBooleanValue()) {
    			command.add("-withzyg");
    		}
    	}
    	// convert2annovar.pl -format pileup -outfile variant.query variant.pileup
    	
    	command.add(outfile);
    	command.add(path2inFile);
    	
    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2outfile),
    			(FileCell) FileCellFactory.create(installPath)};
    	
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
    	m_installpath.saveSettingsTo(settings);
    	m_infile.saveSettingsTo(settings);
    	m_snpqual.saveSettingsTo(settings);
    	m_allallele.saveSettingsTo(settings);
    	m_allelicfrac.saveSettingsTo(settings);
    	m_coverage.saveSettingsTo(settings);
    	m_ifmaxcoverage.saveSettingsTo(settings);
    	m_includeinfo.saveSettingsTo(settings);
    	m_maxcoverage.saveSettingsTo(settings);
    	m_snppvalue.saveSettingsTo(settings);
    	m_species.saveSettingsTo(settings);
    	m_withzyg.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_installpath.loadSettingsFrom(settings);
    	m_infile.loadSettingsFrom(settings);
    	m_snpqual.loadSettingsFrom(settings);
    	m_allallele.loadSettingsFrom(settings);
    	m_allelicfrac.loadSettingsFrom(settings);
    	m_coverage.loadSettingsFrom(settings);
    	m_ifmaxcoverage.loadSettingsFrom(settings);
    	m_includeinfo.loadSettingsFrom(settings);
    	m_maxcoverage.loadSettingsFrom(settings);
    	m_snppvalue.loadSettingsFrom(settings);
    	m_species.loadSettingsFrom(settings);
    	m_withzyg.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_installpath.validateSettings(settings);
    	m_infile.validateSettings(settings);
    	m_snpqual.validateSettings(settings);
    	m_allallele.validateSettings(settings);
    	m_allelicfrac.validateSettings(settings);
    	m_coverage.validateSettings(settings);
    	m_ifmaxcoverage.validateSettings(settings);
    	m_includeinfo.validateSettings(settings);
    	m_maxcoverage.validateSettings(settings);
    	m_snppvalue.validateSettings(settings);
    	m_species.validateSettings(settings);
    	m_withzyg.validateSettings(settings);
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

