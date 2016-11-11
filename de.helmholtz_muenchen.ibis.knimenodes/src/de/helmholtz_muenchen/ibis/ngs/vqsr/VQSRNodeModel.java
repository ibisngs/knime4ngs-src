/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.vqsr;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

/**
 * This is the model implementation of VQSR.
 * 
 *
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class VQSRNodeModel extends HTExecutorNodeModel {
    
	
	/**
	 * Logger
	 */
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(VQSRNodeModel.class);
	
	/**
	 * CFGKEYS
	 */
	
    public static final String CFGKEY_GATK						= "gatk";
    public static final String CFGKEY_REF_GENOME				= "refgenome";
    public static final String CFGKEY_XMX						= "xmx";
    public static final String CFGKEY_MODE						= "MODE";
    public static final String CFGKEY_AN						= "AN";
    public static final String CFGKEY_TRANCHE					= "Tranche";
    public static final String CFGKEY_GAUSS						= "Gaussians";
    public static final String CFGKEY_NT						= "NT";

    public static final String CFGKEY_RESOURCE_HAPMAP			= "resource_hapmap";
    public static final String CFGKEY_RESOURCE_OMNI				= "resource_omni";
    public static final String CFGKEY_RESOURCE_1000G			= "resource_1000G";
    public static final String CFGKEY_RESOURCE_DBSNP			= "resource_dbSNP";
    public static final String CFGKEY_RESOURCE_MILLS			= "resource_mills";

    public static final String CFGKEY_RESOURCES_BOOLEAN_HAPMAP	= "resources_BOOLEAN_HAPMAP";
    public static final String CFGKEY_RESOURCES_BOOLEAN_OMNI	= "resources_BOOLEAN_OMNI";
    public static final String CFGKEY_RESOURCES_BOOLEAN_1000G	= "resources_BOOLEAN_1000G";
    public static final String CFGKEY_RESOURCES_BOOLEAN_DBSNP	= "resources_BOOLEAN_DBSNP";
    public static final String CFGKEY_RESOURCES_BOOLEAN_MILLS	= "resources_BOOLEAN_MILLS";

    public static final String CFGKEY_RESOURCES_FILE_HAPMAP			= "resources_FILE_HAPMAP";
    public static final String CFGKEY_RESOURCES_FILE_OMNI			= "resources_FILE_OMNI";
    public static final String CFGKEY_RESOURCES_FILE_1000G			= "resources_FILE_1000G";
    public static final String CFGKEY_RESOURCES_FILE_DBSNP			= "resources_FILE_DBSNP";
    public static final String CFGKEY_RESOURCES_FILE_MILLS			= "resources_FILE_5";

    public static final String CFGKEY_TS_FILTER					= "TS_FILTER";
    
    public static final String CFGKEY_OPT_VAR_RECAL 			= "var_recal_opts";
    public static final String CFGKEY_OPT_APPLY_RECAL			= "apply_recal_opts";
    
    public static final String DEFAULT_MODE						= "SNP";
    public static final String DEFAULT_SNP_AN					= "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -an MQ";
    public static final String DEFAULT_INDEL_AN					= "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff";
    public static final String DEFAULT_TRANCHES					= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0";
    public static final double DEFAULT_TS_FILTER_SNP			= 99.5;
    public static final double DEFAULT_TS_FILTER_INDEL			= 99.0;
    
	/**
	 * Models
	 */
    
    private final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_XMX = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_XMX, 8, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, VQSRNodeModel.DEFAULT_MODE);
    private final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, VQSRNodeModel.DEFAULT_TRANCHES);
    private final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, VQSRNodeModel.DEFAULT_SNP_AN);
    private final SettingsModelIntegerBounded m_GAUSSIANS = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_GAUSS,8,1,Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

    private final SettingsModelString m_RESOURCES_STRING_HAPMAP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_HAPMAP, "hapmap,known=false,training=true,truth=true,prior=15.0");
    private final SettingsModelString m_RESOURCES_STRING_OMNI = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_OMNI, "omni,known=false,training=true,truth=true,prior=12.0");
    private final SettingsModelString m_RESOURCES_STRING_1000G = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_1000G, "1000G,known=false,training=true,truth=false,prior=10.0");
    private final SettingsModelString m_RESOURCES_STRING_DBSNP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_DBSNP, "dbsnp,known=true,training=false,truth=false,prior=2.0");
    private final SettingsModelString m_RESOURCES_STRING_MILLS = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_MILLS, "mills,known=false,training=true,truth=true,prior=12.0");
    
    //TODO remove for publication
    private final SettingsModelString m_RESOURCES_FILE_HAPMAP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_HAPMAP,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/hapmap_3.3.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_OMNI = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_OMNI,"");    // "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_omni2.5.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_1000G = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1000G,"");  // "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_phase1.snps.high_confidence.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_DBSNP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_DBSNP,"");  // "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/dbsnp_138.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_MILLS = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_MILLS,"");  // "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf");

    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_HAPMAP = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_HAPMAP, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_OMNI = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_OMNI, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_1000G = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1000G, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_DBSNP = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_DBSNP, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_MILLS = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_MILLS, false);

    private final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,VQSRNodeModel.DEFAULT_TS_FILTER_SNP,1,100);

    private final SettingsModelOptionalString m_OPT_VAR_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_VAR_RECAL,"",false);
    private final SettingsModelOptionalString m_OPT_APPLY_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_APPLY_RECAL,"",false);
    
	public static final String OUT_COL1 = "VQSR VARIANTS";
	private int vcf_index;
    
    /**
     * Constructor for the node model.
     */
    protected VQSRNodeModel() {
    
        super(1, 1);
        
        addSetting(m_AN);
        addSetting(m_GATK);
        addSetting(m_MODE);
        addSetting(m_NT);
        addSetting(m_GAUSSIANS);
        addSetting(m_REF_GENOME);
        addSetting(m_TRANCHE);
        addSetting(m_XMX);
        
        addSetting(m_RESOURCES_BOOLEAN_HAPMAP);
        addSetting(m_RESOURCES_BOOLEAN_OMNI);
        addSetting(m_RESOURCES_BOOLEAN_1000G);
        addSetting(m_RESOURCES_BOOLEAN_DBSNP);
        addSetting(m_RESOURCES_BOOLEAN_MILLS);
        addSetting(m_RESOURCES_FILE_HAPMAP);
        addSetting(m_RESOURCES_FILE_OMNI);
        addSetting(m_RESOURCES_FILE_1000G);
        addSetting(m_RESOURCES_FILE_DBSNP);
        addSetting(m_RESOURCES_FILE_MILLS);
        addSetting(m_RESOURCES_STRING_HAPMAP);
        addSetting(m_RESOURCES_STRING_OMNI);
        addSetting(m_RESOURCES_STRING_1000G);
        addSetting(m_RESOURCES_STRING_DBSNP);
        addSetting(m_RESOURCES_STRING_MILLS);
        addSetting(m_TS_FILTER);
        
        addSetting(m_OPT_VAR_RECAL);
        addSetting(m_OPT_APPLY_RECAL);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String PATH2GATK = m_GATK.getStringValue();
    	String PATH2REFSEQ = m_REF_GENOME.getStringValue();
    	
    	String INFILE = inData[0].iterator().next().getCell(vcf_index).toString();
    	if(CompatibilityChecker.inputFileNotOk(INFILE)) {
    		throw new InvalidSettingsException("No VCF file in the input table or VCF file does not exist!");
    	}
    	
    	String outFile = IO.replaceFileExtension(INFILE, m_MODE.getStringValue()+"_VQSR.vcf");
    	
    	File lockFile1 = new File(IO.replaceFileExtension(outFile,"var_recal" + SuccessfulRunChecker.LOCK_ENDING));
    	File lockFile2 = new File(IO.replaceFileExtension(outFile,"apply_recal" + SuccessfulRunChecker.LOCK_ENDING));
    	String stdOut1 = IO.replaceFileExtension(outFile,"var_recal.stdOut");
    	String stdErr1 = IO.replaceFileExtension(outFile,"var_recal.stdErr");
    	String stdOut2 = IO.replaceFileExtension(outFile,"apply_recal.stdOut");
    	String stdErr2 = IO.replaceFileExtension(outFile,"apply_recal.stdErr");
    	
    	super.executeCommand(createRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,lockFile1,stdOut1,stdErr1);
    	super.executeCommand(applyRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,lockFile2,stdOut2,stdErr2);

    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(outFile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     *Creates Recalibration Command
     * @return
     */
    private String[] createRecalibrationCommand(String PATH2GATK, String PATH2REFSEQ, String INFILE){
    	
    	ArrayList<String> command = new ArrayList<String>();
    
    	String recalFile 	= IO.replaceFileExtension(INFILE, m_MODE.getStringValue()+ "_VQSR.recal");
    	String tranchesFile	= IO.replaceFileExtension(INFILE, m_MODE.getStringValue() +"_VQSR.tranches");
    	String plotFile		= IO.replaceFileExtension(INFILE, m_MODE.getStringValue()+"_VQSR.plots.R");
    	
    	int memory = m_XMX.getIntValue() * m_NT.getIntValue();
    	
    	command.add("java -Xmx"+memory+"G -jar");
    	command.add(PATH2GATK);
    	command.add("-T VariantRecalibrator");
    	command.add("-R "+PATH2REFSEQ);
    	command.add("-input "+INFILE);
    	
    	//Add resources
    	if(m_RESOURCES_BOOLEAN_HAPMAP.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_HAPMAP.getStringValue()+" "+m_RESOURCES_FILE_HAPMAP.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_OMNI.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_OMNI.getStringValue()+" "+m_RESOURCES_FILE_OMNI.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_1000G.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_1000G.getStringValue()+" "+m_RESOURCES_FILE_1000G.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_DBSNP.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_DBSNP.getStringValue()+" "+m_RESOURCES_FILE_DBSNP.getStringValue());
    	}
    	
    	if(m_RESOURCES_BOOLEAN_MILLS.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_MILLS.getStringValue()+" "+m_RESOURCES_FILE_MILLS.getStringValue());
    	}
    	
    	command.add(m_AN.getStringValue());
    	
    	command.add("-mode "+m_MODE.getStringValue());
    	
    	command.add("-recalFile "+recalFile);
    	command.add("-tranchesFile "+tranchesFile);
    	command.add("-rscriptFile "+plotFile);
    	command.add("-nt "+m_NT.getIntValue());
    	command.add(m_TRANCHE.getStringValue());
    	command.add("--maxGaussians "+m_GAUSSIANS.getIntValue());
    	
    	if(m_OPT_VAR_RECAL.isActive()) {
    		command.add(m_OPT_VAR_RECAL.getStringValue());
    	}
    	
    	return new String[]{StringUtils.join(command, " ")};
    }
    
    /**
     *Creates Recalibration Command
     * @return
     */
    private String[] applyRecalibrationCommand(String PATH2GATK, String PATH2REFSEQ,String INFILE){
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	//Process Infile Name to obtain output file names
    	String infile 		= INFILE;
    	String recalFile 	= IO.replaceFileExtension(INFILE, m_MODE.getStringValue()+ "_VQSR.recal");
    	String tranchesFile	= IO.replaceFileExtension(INFILE, m_MODE.getStringValue() +"_VQSR.tranches");
    	String outFile = IO.replaceFileExtension(INFILE, m_MODE.getStringValue()+"_VQSR.vcf");

    	command.add("java -Xmx"+m_XMX.getIntValue()+"G  -jar");
    	command.add(PATH2GATK);
    	command.add("-T ApplyRecalibration");
    	command.add("-R "+PATH2REFSEQ);
    	command.add("-input "+infile); 	
       	command.add("-mode "+m_MODE.getStringValue());
    	
    	command.add("-recalFile "+recalFile);
    	command.add("-tranchesFile "+tranchesFile);
    	command.add("-o "+outFile);

    	command.add("-ts_filter_level "+m_TS_FILTER.getDoubleValue());
    	if(m_OPT_APPLY_RECAL.isActive()) {
    		command.add(m_OPT_APPLY_RECAL.getStringValue());
    	}
    	    	
    	return new String[]{StringUtils.join(command, " ")};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_GATK.getStringValue())) {
    		throw new InvalidSettingsException("Set path to the GenomeAnalysisTK.jar!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_REF_GENOME.getStringValue())) {
    		throw new InvalidSettingsException("Set reference genome!");
    	}
    	
    	if(m_RESOURCES_BOOLEAN_HAPMAP.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_RESOURCES_FILE_HAPMAP.getStringValue())){
    		throw new InvalidSettingsException("Set HapMap reference data set!");
    	}
    	if(m_RESOURCES_BOOLEAN_OMNI.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_RESOURCES_FILE_OMNI.getStringValue())){
    		throw new InvalidSettingsException("Set Omni reference data set!");
    	}
    	if(m_RESOURCES_BOOLEAN_1000G.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_RESOURCES_FILE_1000G.getStringValue())){
    		throw new InvalidSettingsException("Set 1000G reference data set!");
    	}
    	if(m_RESOURCES_BOOLEAN_DBSNP.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_RESOURCES_FILE_DBSNP.getStringValue())){
    		throw new InvalidSettingsException("Set dbSNP reference data set!");
    	}
    	if(m_RESOURCES_BOOLEAN_MILLS.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_RESOURCES_FILE_MILLS.getStringValue())){
    		throw new InvalidSettingsException("Set Mills reference data set!");
    	}

    	if(m_AN.getStringValue().equals("")) {
    		throw new InvalidSettingsException("Annotation field must not be empty!");
    	}
    	
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }
}

