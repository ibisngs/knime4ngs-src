package de.helmholtz_muenchen.ibis.ngs.vcfnormalizer;

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

import de.helmholtz_muenchen.ibis.ngs.vep.VEPNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of VCFNormalizer.
 * 
 *
 * @author tim.jeske
 */
public class VCFNormalizerNodeModel extends NodeModel {
    
	// the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(VEPNodeModel.class);
    
    public static boolean optionalPort=false;
    
    static final String CFGKEY_VT = "vt_path";
    final SettingsModelString m_vt_path = new SettingsModelString(CFGKEY_VT,"");
    
    static final String CFGKEY_VCF = "vcf_infile";
    final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCF,"");
    
    static final String CFGKEY_REF_GENOME = "ref_genome";
    final SettingsModelString m_ref_genome = new SettingsModelString(CFGKEY_REF_GENOME,"");
    
    public static final String OUT_COL = "Path2normalizedVCF";
	
    /**
     * Constructor for the node model.
     */
    protected VCFNormalizerNodeModel() {
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String vt = m_vt_path.getStringValue();
    	String vcfin = m_vcfin.getStringValue();
    	String ref_genome = m_ref_genome.getStringValue();
    	
    	String outfile = vcfin.replace(".vcf.gz", ".vcf");
    	outfile = outfile.replace(".vcf", ".normalized.vcf");
    	
    	String cmd="";
    	
    	if(vcfin.endsWith(".vcf.gz")) {
    		cmd = "zless "+vcfin;
    	} else if (vcfin.endsWith(".vcf")) {
    		cmd = "less "+vcfin;
    	} else {
    		logger.error("VCF infile not correctly specified!");
    	}
    	
    	cmd += " | sed 's/ID=AD,Number=./ID=AD,Number=R/'";
    	
    	if(vt.equals("") || Files.notExists(Paths.get(vt))) {
    		logger.error("vt directory not correctly specified!");
    	} else {
    		if(vt.endsWith("/")) {
    			vt = vt+ "vt";
    		} else {
    			vt = vt + "/vt";
    		}
    		cmd += " | "+ vt + " decompose -s -"; 
    	}
    	
    	if(ref_genome.equals("") || Files.notExists(Paths.get(ref_genome))) {
    		logger.error("Refgenome not correctly specified!");
    	} else {
    		cmd += " | " + vt + " normalize -r "+ref_genome+" - > "+ outfile;
    	}
    	
    	Executor.executeCommand(new String[]{"/usr/bin/bash","-c",cmd}, exec, logger);
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile)};
    	
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
    					new DataColumnSpecCreator(OUT_COL, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_vt_path.saveSettingsTo(settings);
         m_vcfin.saveSettingsTo(settings);
         m_ref_genome.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_vt_path.loadSettingsFrom(settings);
        m_vcfin.loadSettingsFrom(settings);
        m_ref_genome.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_vt_path.validateSettings(settings);
        m_vcfin.validateSettings(settings);
        m_ref_genome.validateSettings(settings);
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

