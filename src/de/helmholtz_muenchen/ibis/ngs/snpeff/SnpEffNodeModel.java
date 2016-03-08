package de.helmholtz_muenchen.ibis.ngs.snpeff;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of SnpEff.
 * 
 *
 * @author Sebastian Kopetzky
 */
public class SnpEffNodeModel extends HTExecutorNodeModel {
    /*Mandatory options*/
	
	public static final String CFGKEY_SNPEFF_FOLDER = "snpeff_folder";
	public static final String CFGKEY_VCF_FILE = "vcf_file";
	public static final String CFGKEY_DATABASE = "database";

	//directory containing the snpEff.jar and scripts
	private final SettingsModelString m_snpeff_folder = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_SNPEFF_FOLDER,"");
	//name of the database to use
	private final SettingsModelString m_database = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_DATABASE,"");
	//input vcf file
	private final SettingsModelString m_vcf_file = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_VCF_FILE,"");

	/*Sequence change filter options*/
	public static final String CFGKEY_USEMINQ = "use_minq";
	public static final String CFGKEY_MINQ = "minq";
	public static final String CFGKEY_USEMINC = "use_minc";
	public static final String CFGKEY_MINC = "minc";
	public static final String CFGKEY_DEL = "del";
	public static final String CFGKEY_INS = "ins";
	public static final String CFGKEY_HOM = "hom";
	public static final String CFGKEY_HET = "het";
	public static final String CFGKEY_MNP = "mnp";
	public static final String CFGKEY_SNP = "snp";
	
	//Filter out variants with quality lower minQ	
	private final SettingsModelBoolean m_use_minq = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEMINQ, false);
	private final SettingsModelDoubleBounded m_minq = new SettingsModelDoubleBounded(
			SnpEffNodeModel.CFGKEY_MINQ, 0.0, 0.0, Double.MAX_VALUE);
	//Filter out variants with coverage lower minC
	private final SettingsModelBoolean m_use_minc = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEMINC, false);
	private final SettingsModelIntegerBounded m_minc = new SettingsModelIntegerBounded(
			SnpEffNodeModel.CFGKEY_MINC, 0, 0, Integer.MAX_VALUE);
	//Analyze deletions only
	private final SettingsModelBoolean m_del = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_DEL, false);
	//Analyze insertions only
	private final SettingsModelBoolean m_ins = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_INS, false);
	//Analyze homozygous variants only
	private final SettingsModelBoolean m_hom = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_HOM, false);
	//Analyze heterozygous variants only
	private final SettingsModelBoolean m_het = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_HET, false);
	//Analyze only MNPs (multiple nucleotide polymorphisms)
	private final SettingsModelBoolean m_mnp = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_MNP, false);
	//Only SNPs (single nucleotide polymorphisms)
	private final SettingsModelBoolean m_snp = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_SNP, false);

	
	/*Results filter options*/
	public static final String CFGKEY_USEBEDFILE="usebedfile";
	public static final String CFGKEY_BED_FILE = "bed_file";
	
	public static final String CFGKEY_NO_DOWNSTREAM = "no_downstream";
	public static final String CFGKEY_NO_INTERGENIC = "no_intergenic";
	public static final String CFGKEY_NO_INTRONIC = "no_intronic";
	public static final String CFGKEY_NO_UPSTREAM = "no_upstream";
	public static final String CFGKEY_NO_UTR = "no_utr";

	//checkbox if bed filter should be used
	private final SettingsModelBoolean m_usebedfile = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEBEDFILE, false);
	//bed file path
	private final SettingsModelString m_bed_file = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_BED_FILE,"");
	
	private final SettingsModelBoolean m_no_downstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_DOWNSTREAM, false);
	private final SettingsModelBoolean m_no_intergenic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTERGENIC, false);
	private final SettingsModelBoolean m_no_intronic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTRONIC, false);
	private final SettingsModelBoolean m_no_upstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UPSTREAM, false);
	private final SettingsModelBoolean m_no_utr = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UTR, false);
	
	/*Annotations options*/
	public static final String CFGKEY_LOF = "lof";
	
	private final SettingsModelBoolean m_lof = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_LOF, false);
	
	/*Other variables*/
	private boolean optionalPort = false;
	
	//The Output Col Names
	public static final String OUT_COL1 = "snpEffDirectory";
	public static final String OUT_COL2 = "outputFile";
	
	/**
     * Constructor for the node model.
     */
    protected SnpEffNodeModel() {
    
        super(OptionalPorts.createOPOs(1,1), OptionalPorts.createOPOs(1));
        
        /*Mandatory options*/
    	addSetting(m_snpeff_folder);
    	addSetting(m_database);
    	addSetting(m_vcf_file);
    	
    	/*Sequence change filter options*/
    	addSetting(m_use_minq);
    	addSetting(m_minq);
    	addSetting(m_use_minc);
    	addSetting(m_minc);
    	
    	addSetting(m_del);
    	addSetting(m_ins);
    	addSetting(m_hom);
    	addSetting(m_het);
    	addSetting(m_mnp);
    	addSetting(m_snp);
    	
    	/*Results filter options*/
    	addSetting(m_usebedfile);
    	addSetting(m_bed_file);
    	
    	addSetting(m_no_downstream);
    	addSetting(m_no_intergenic);
    	addSetting(m_no_intronic);
    	addSetting(m_no_upstream);
    	addSetting(m_no_utr);
    	
    	/*Annotations options*/
    	addSetting(m_lof);
        
        m_bed_file.setEnabled(false);
        
        m_minc.setEnabled(false);
        m_minq.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	/*Get parameters*/
    	String snpEffDirectory;
    	String dbName;
    	String vcfFile = m_vcf_file.getStringValue();
    	//TODO: add parameter for outFile?
    	String outFile = vcfFile.substring(0, vcfFile.length()-4);
    	
    	if(optionalPort){
    		//get param from port
    		snpEffDirectory = inData[0].iterator().next().getCell(0).toString();
    		dbName = inData[0].iterator().next().getCell(1).toString();
    	}
    	else{
    		//get param normally
    		snpEffDirectory = m_snpeff_folder.getStringValue();
    		dbName = m_database.getStringValue();
    	}
    	   	
    	/*Make & run command*/
    	//TODO change -Xmx value? Dynamically?
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-Xmx2500m -jar "+snpEffDirectory+"/snpEff.jar eff");
    	command.add("-s "+outFile + "_summary.html");
    	command.add("-v "+dbName);
    	
    	//Sequence change filter options
    	if(m_use_minq.getBooleanValue()){
    		command.add("-minQ " + new Double(m_minq.getDoubleValue()).toString());
    	}
    	if(m_use_minc.getBooleanValue()){
    		command.add("-minC " + new Integer(m_minc.getIntValue()).toString());
    	}
    	if(m_del.getBooleanValue()){
    		command.add("-del");
    	}
    	else if(m_ins.getBooleanValue()){
    		command.add("-ins");
    	}
    	if(m_het.getBooleanValue()){
    		command.add("-het");
    	}
    	else if(m_hom.getBooleanValue()){
    		command.add("-hom");
    	}
    	if(m_snp.getBooleanValue()){
    		command.add("-snp");
    	}
    	else if(m_mnp.getBooleanValue()){
    		command.add("-mnp");
    	}
    	
    	//Result filter options
    	if(m_usebedfile.getBooleanValue()){
    		command.add("-fi " + m_bed_file.getStringValue());
    	}
    	if(m_no_downstream.getBooleanValue()){
    		command.add("-no-downstream");
    	}
    	if(m_no_intergenic.getBooleanValue()){
    		command.add("-no-intergenic");
    	}
    	if(m_no_intronic.getBooleanValue()){
    		command.add("-no-intron");
    	}
    	if(m_no_upstream.getBooleanValue()){
    		command.add("-no-upstream");
    	}
    	if(m_no_utr.getBooleanValue()){
    		command.add("-no-utr");
    	}
    	
    	/*Annotations options*/
    	if(m_lof.getBooleanValue()){
    		command.add("-lof");
    	}
    	
    	command.add(vcfFile);

    	
    	/**Execute**/
    	String lockFile = outFile + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile),outFile + ".eff.vcf");
    	
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,outFile + ".eff.vcf");
//    	logBuffer.append(ShowOutput.getNodeEndTime());
//    	ShowOutput.writeLogFile(logBuffer);
    	
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(snpEffDirectory),
    			(FileCell) FileCellFactory.create(outFile + ".eff.vcf")};
    	
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

		//Check OptionalInputPort
		try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			String[] colnames = inSpecs[0].getColumnNames();
			if(colnames[0].equals("snpEffDirectory") && colnames[1].equals("database")){
				m_snpeff_folder.setEnabled(false);
				m_database.setEnabled(false);
			}
			else{
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			}

			//m_bfast.setEnabled(false);
    	}catch(NullPointerException npe){
			m_snpeff_folder.setEnabled(true);
			m_database.setEnabled(true);
    	}
		
        if(!optionalPort && m_snpeff_folder.getStringValue().length() == 0){
        	throw new InvalidSettingsException("Specify path to snpEff!");
        }
        
        if(!optionalPort && m_database.getStringValue().length() == 0){
        	throw new InvalidSettingsException("Specify the database to use!");
        }
        
        if(m_vcf_file.getStringValue().length() == 0){
        	throw new InvalidSettingsException("Specify an input vcf file!");
        }
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
    }

//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void saveSettingsTo(final NodeSettingsWO settings) {
//    	
//    	super.saveSettingsTo(settings);
//    	
//    	/*Mandatory options*/
//    	m_snpeff_folder.saveSettingsTo(settings);
//    	m_database.saveSettingsTo(settings);
//    	m_vcf_file.saveSettingsTo(settings);
//    	
//    	/*Sequence change filter options*/
//    	m_use_minq.saveSettingsTo(settings);
//    	m_minq.saveSettingsTo(settings);
//    	m_use_minc.saveSettingsTo(settings);
//    	m_minc.saveSettingsTo(settings);
//    	
//    	m_del.saveSettingsTo(settings);
//    	m_ins.saveSettingsTo(settings);
//    	m_hom.saveSettingsTo(settings);
//    	m_het.saveSettingsTo(settings);
//    	m_mnp.saveSettingsTo(settings);
//    	m_snp.saveSettingsTo(settings);
//    	
//    	/*Results filter options*/
//    	m_usebedfile.saveSettingsTo(settings);
//    	m_bed_file.saveSettingsTo(settings);
//    	
//    	m_no_downstream.saveSettingsTo(settings);
//    	m_no_intergenic.saveSettingsTo(settings);
//    	m_no_intronic.saveSettingsTo(settings);
//    	m_no_upstream.saveSettingsTo(settings);
//    	m_no_utr.saveSettingsTo(settings);
//    	
//    	/*Annotations options*/
//    	m_lof.saveSettingsTo(settings);
//
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//        
//    	super.loadValidatedSettingsFrom(settings);
//    	
//    	/*Mandatory options*/
//    	m_snpeff_folder.loadSettingsFrom(settings);
//    	m_database.loadSettingsFrom(settings);
//    	m_vcf_file.loadSettingsFrom(settings);
//    	
//    	/*Sequence change filter options*/
//    	m_use_minq.loadSettingsFrom(settings);
//    	m_minq.loadSettingsFrom(settings);
//    	m_use_minc.loadSettingsFrom(settings);
//    	m_minc.loadSettingsFrom(settings);
//    	m_mnp.loadSettingsFrom(settings);
//    	m_snp.loadSettingsFrom(settings);
//    	
//    	m_del.loadSettingsFrom(settings);
//    	m_ins.loadSettingsFrom(settings);
//    	m_hom.loadSettingsFrom(settings);
//    	m_het.loadSettingsFrom(settings);
//    	
//    	/*Results filter options*/
//    	m_usebedfile.loadSettingsFrom(settings);
//    	m_bed_file.loadSettingsFrom(settings);
//    	m_no_downstream.loadSettingsFrom(settings);
//    	m_no_intergenic.loadSettingsFrom(settings);
//    	m_no_intronic.loadSettingsFrom(settings);
//    	m_no_upstream.loadSettingsFrom(settings);
//    	m_no_utr.loadSettingsFrom(settings);
//    	
//    	/*Annotations options*/
//    	m_lof.loadSettingsFrom(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void validateSettings(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//
//    	super.validateSettings(settings);
//    	
//    	/*Mandatory options*/
//    	m_snpeff_folder.validateSettings(settings);
//    	m_database.validateSettings(settings);
//    	m_vcf_file.validateSettings(settings);
//    	
//    	/*Sequence change filter options*/
//    	m_use_minq.validateSettings(settings);
//    	m_minq.validateSettings(settings);
//    	m_use_minc.validateSettings(settings);
//    	m_minc.validateSettings(settings);
//    	
//    	m_del.validateSettings(settings);
//    	m_ins.validateSettings(settings);
//    	m_hom.validateSettings(settings);
//    	m_het.validateSettings(settings);
//    	m_mnp.validateSettings(settings);
//    	m_snp.validateSettings(settings);
//    	
//    	/*Results filter options*/
//    	m_usebedfile.validateSettings(settings);
//    	m_bed_file.validateSettings(settings);
//    	m_no_downstream.validateSettings(settings);
//    	m_no_intergenic.validateSettings(settings);
//    	m_no_intronic.validateSettings(settings);
//    	m_no_upstream.validateSettings(settings);
//    	m_no_utr.validateSettings(settings);
//    	
//    	/*Annotations options*/
//    	m_lof.validateSettings(settings);
//    	
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

