package de.helmholtz_muenchen.ibis.ngs.reordervcf;

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

import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


/**
 * This is the model implementation of ReorderVCF.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class ReorderVCFNodeModel extends HTExecutorNodeModel {
    
	public static final String CFGKEY_REFERENCE_VCF = "REFERENCE_VCF";
	public static final String CFGKEY_VCFTOOLS 		= "VCFTOOLS";
	
	private final SettingsModelString m_REFERENCE_VCF 		= new SettingsModelString(ReorderVCFNodeModel.CFGKEY_REFERENCE_VCF, "---");
	private final SettingsModelString m_VCFTOOLS 			= new SettingsModelString(ReorderVCFNodeModel.CFGKEY_VCFTOOLS, "---");
	
	//The Output Col Names
	public static final String OUT_COL1 = "REORDEREDVCF";
	
    /**
     * Constructor for the node model.
     */
    protected ReorderVCFNodeModel() {
        super(1, 1);
        addSetting(m_REFERENCE_VCF);
        addSetting(m_VCFTOOLS);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {


    	String VCFToOrder 	= inData[0].iterator().next().getCell(0).toString();
    	String ReferenceVCF	= m_REFERENCE_VCF.getStringValue();
    	String VCFTOOLS		= m_VCFTOOLS.getStringValue();
    	String REORDERED	= VCFToOrder.replace(".vcf", ".reordered.vcf");
    	
    	boolean bgzipINFILE = false;
    	boolean tabixINFILE = false;
    	boolean bgzipREFVCF = false;
    	boolean tabixREFVCF = false;
    	
    	//Check Infiles and determine what has to be done
    	if(!VCFToOrder.endsWith(".vcf.gz")){bgzipINFILE=true;}
    	if(!new File(VCFToOrder+".vcf.gz.tbi").exists()){tabixINFILE=true;}
    	if(!ReferenceVCF.endsWith(".vcf.gz")){bgzipREFVCF=true;}
    	if(!new File(ReferenceVCF+".vcf.gz.tbi").exists()){tabixREFVCF=true;}    	
    	
    	//Prepare Infile
    	if(bgzipINFILE){
    		VCFToOrder = bgzip(VCFToOrder, exec);
    		if(tabixINFILE){
    			tabix(VCFToOrder, exec);
    		}
    	}else{
    		if(tabixINFILE){
    			tabix(VCFToOrder, exec);
    		}
    	}
    	//Prepare Reference VCF File
    	if(bgzipREFVCF){
    		ReferenceVCF = bgzip(ReferenceVCF, exec);
    		if(tabixREFVCF){
    			tabix(ReferenceVCF, exec);
    		}
    	}else{
    		if(tabixREFVCF){
    			tabix(ReferenceVCF, exec);
    		}
    	}    	
    	
    	//REORDER
     	ArrayList<String> command = new ArrayList<String>();
    	command = new ArrayList<String>();
    	command.add("perl");
    	command.add(VCFTOOLS+"/vcf-shuffle-cols");
    	command.add("-t "+ReferenceVCF);
    	command.add(VCFToOrder);
    	
    	File lockFile = new File(REORDERED + SuccessfulRunChecker.LOCK_ENDING);
    	
    	super.executeCommand(new String[]{StringUtils.join(command, " ")},exec,new String[]{"PERL5LIB=PERL5LIB:"+VCFTOOLS},lockFile,REORDERED,REORDERED+".stdErr",null,null,null);
    	

    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(REORDERED)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();

        return new BufferedDataTable[]{outTable};
    }

    /**
     * Execute bgzip on Infile
     * @param INFILE
     * @param exec
     * @return Outfile
     * @throws Exception 
     */
    private String bgzip(String INFILE,final ExecutionContext exec) throws Exception{
    	String GZFILE = INFILE+".gz";
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("/home/software/bin/bgzip -c "+INFILE);
    	
    	File lockFile = new File(GZFILE + SuccessfulRunChecker.LOCK_ENDING);
    	
    	super.executeCommand(new String[]{StringUtils.join(command, " ")},exec,lockFile,GZFILE);
    	
    	return GZFILE;
    }
    
    /**
     * Execute tabix on Infile
     * @param INFILE
     * @param exec
     * @throws Exception 
     */
    private void tabix(String INFILE,final ExecutionContext exec) throws Exception{
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("/home/software/bin/tabix -p vcf "+INFILE);
    	File lockFile = new File(INFILE+".tbi" + SuccessfulRunChecker.LOCK_ENDING);
    	super.executeCommand(new String[]{StringUtils.join(command, " ")},exec,lockFile);
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
        
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void saveSettingsTo(final NodeSettingsWO settings) {
//    	/** added for HTE **/
//    	super.saveSettingsTo(settings);
//    	
//         m_REFERENCE_VCF.saveSettingsTo(settings);
//         m_VCFTOOLS.saveSettingsTo(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//    	/** added for HTE **/
//    	super.loadValidatedSettingsFrom(settings);
//    	
//        m_REFERENCE_VCF.loadSettingsFrom(settings);
//        m_VCFTOOLS.loadSettingsFrom(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void validateSettings(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//    	/** added for HTE **/
//    	super.validateSettings(settings);
//    	
//        m_REFERENCE_VCF.validateSettings(settings);
//        m_VCFTOOLS.validateSettings(settings);
//    }
    
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

