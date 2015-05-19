package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of LOFSummarizer.
 * 
 *
 * @author tim.jeske
 */
public class LOFSummaryNodeModel extends NodeModel {
    
	// the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(LOFSummaryNodeModel.class);
	
	//input file
	static final String CFGKEY_VCF_INFILE = "vcf_infile";
	final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCF_INFILE,"");
	
	static final String CFGKEY_CDS_INFILE = "cds_infile";
	final SettingsModelString m_cdsin = new SettingsModelString(CFGKEY_CDS_INFILE,"");
	
	//selected annotation
    static final String CFGKEY_ANNOTATION="annotation";
    static final String[] ANNOTATIONS_AVAILABLE={"VAT","VEP"};
    final SettingsModelString m_annotation = new SettingsModelString(CFGKEY_ANNOTATION, "");
	
	//output
	public static final String OUT_COL1 = "Path2LOF_Statistics";
	public static final String OUT_COL2 = "Path2Gene_Statistics";
	public static final String OUT_COL3 = "Path2Sample_Statistics";
	public static final String OUT_COL4 = "Path2Transcript_Statistics";
	
	public static boolean optionalPort=false;
	
    /**
     * Constructor for the node model.
     */
    protected LOFSummaryNodeModel() {
    
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String vcf_infile;
    	
    	if(optionalPort){	//Input Table available
    		//Get File via Table
    		vcf_infile = inData[0].iterator().next().getCell(0).toString();
    	}else{
    		//Get File via FileSelector
    		vcf_infile = m_vcfin.getStringValue();
    		if(vcf_infile.equals("") || Files.notExists(Paths.get(vcf_infile))) {
    			logger.error("No input vcf file specified!");
    		}
    	}
    	logger.info("Infile: "+vcf_infile);

    	String cds_file = m_cdsin.getStringValue();
    	if(cds_file.equals("") || Files.notExists(Paths.get(cds_file))) {
    		logger.error("No CDS file specified! Variant effect (full or partial) cannot be calculated.");
    	}
    	
    	Summarizer summy = null;
    	String annotation = m_annotation.getStringValue();
    	if(annotation.equals("VAT")) {
    		summy = new VATSummarizer(vcf_infile, cds_file);
    	} else if(annotation.equals("VEP")) {
    		summy = new VEPSummarizer(vcf_infile, cds_file);
    	}
    	
    	String LOF_Summary[] = new String[3];
    	LOF_Summary[0] = summy.getLoFStatistic();
    	LOF_Summary[1] = summy.getGeneStatistic();
    	LOF_Summary[2] = summy.getSampleStatistic();
//    	LOF_Summary[3] = summy.getTranscriptStatistic();
    	//TODO add to output table
//    	summy.getAnnotationStatistic();
//    	summy.getKnockOutGenes();
//    	summy.getSampleKOGeneMatrix();
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(LOF_Summary[0]),
    			(FileCell) FileCellFactory.create(LOF_Summary[1]),
    			(FileCell) FileCellFactory.create(LOF_Summary[2])};
    	
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
        optionalPort = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
		}catch(NullPointerException e){}

        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_vcfin.saveSettingsTo(settings);
    	m_cdsin.saveSettingsTo(settings);
    	m_annotation.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_vcfin.loadSettingsFrom(settings);
    	m_cdsin.loadSettingsFrom(settings);
    	m_annotation.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_vcfin.validateSettings(settings);
    	m_cdsin.validateSettings(settings);
    	m_annotation.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}