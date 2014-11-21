package de.helmholtz_muenchen.ibis.ngs.summarizelof;

import java.io.File;
import java.io.IOException;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of SummarizeLOF.
 * 
 *
 * @author 
 */
public class SummarizeLOFNodeModel extends NodeModel {
    
	public static final String CFGKEY_PED_FILE	 	= "PED_File";
	public static final String CFGKEY_VCF_INFILE 	= "VCF_INFILE";
	
	
	private final SettingsModelString m_pedfile = new SettingsModelString(CFGKEY_PED_FILE,"");
	private final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCF_INFILE,"");
	
	//The Output Col Names
	public static final String OUT_COL1 = "LOF SAMPLE SUMMARY";
	public static final String OUT_COL2 = "LOF GENE SUMMARY";
	
	public static boolean optionalPort=false;
	
    /**
     * Constructor for the node model.
     */
    protected SummarizeLOFNodeModel() {
    
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String VCF_INFILE;
    	String PED_FILE;
    	
    	System.out.println(optionalPort);
    	if(optionalPort){	//Input Table available
    		//Get File via Table
    		VCF_INFILE = inData[0].iterator().next().getCell(0).toString();
    	}else{
    		//Get File via FileSelector
    		VCF_INFILE = m_vcfin.getStringValue();
    	}
    	System.out.println("Infile: "+VCF_INFILE);
    	PED_FILE = m_pedfile.getStringValue();
    	
    	//Execute
    	String LOF_Summary[] = new SummarizeLOF().getLOFs(PED_FILE, VCF_INFILE);
    	
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(LOF_Summary[0]),
    			(FileCell) FileCellFactory.create(LOF_Summary[1])};
    	
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
        // TODO: generated method stub
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

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_pedfile.saveSettingsTo(settings);
         m_vcfin.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_pedfile.loadSettingsFrom(settings);
        m_vcfin.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_pedfile.validateSettings(settings);
        m_vcfin.validateSettings(settings);
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

