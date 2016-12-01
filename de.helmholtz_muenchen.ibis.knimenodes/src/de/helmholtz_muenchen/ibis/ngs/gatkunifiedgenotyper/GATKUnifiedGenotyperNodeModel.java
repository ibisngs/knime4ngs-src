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
package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;

import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;


/**
 * This is the model implementation of GATKUnifiedGenotyper.
 * 
 *
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class GATKUnifiedGenotyperNodeModel extends GATKNodeModel {
    

    protected static final NodeLogger logger = NodeLogger.getLogger(GATKUnifiedGenotyperNodeModel.class);
    
    //output folder
    static final String CFGKEY_OUTFOLDER = "outfolder";
    private final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
    // use dbsnp snp set
 	static final String CFGKEY_USE_DBSNP="use_dbsnp";
 	static final boolean DEF_USE_DBSNP=true;
 	private final SettingsModelBoolean m_use_dbsnp=new SettingsModelBoolean(CFGKEY_USE_DBSNP, DEF_USE_DBSNP);
 	// path to dbsnp snp set
 	static final String CFGKEY_DBSNP_FILE="dbsnp_file";
 	static final String DEF_DBSNP_FILE="";
 	private final SettingsModelString m_dbsnp_file=new SettingsModelString(CFGKEY_DBSNP_FILE, DEF_DBSNP_FILE);
 	// call snps, indels or both
 	static final String CFGKEY_VARIANT_TYPE="variant_type";
 	static final String [] AVAIL_VARIANT_TYPE={"SNP", "INDEL", "BOTH"};
 	static final String DEF_VARIANT_TYPE=AVAIL_VARIANT_TYPE[2];
	private final SettingsModelString m_variant_type= new SettingsModelString(CFGKEY_VARIANT_TYPE, DEF_VARIANT_TYPE);
    // number of threads for target creator
    static final String CFGKEY_NUM_THREADS="num_threads";
    static final int DEF_NUM_THREADS=1;
    static final int MIN_NUM_THREADS=1;
    static final int MAX_NUM_THREADS=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_num_threads= new SettingsModelIntegerBounded(CFGKEY_NUM_THREADS, DEF_NUM_THREADS, MIN_NUM_THREADS, MAX_NUM_THREADS);
    
    // additional options
    
    // confidence threshold for calling
    static final String CFGKEY_CALL_MIN_CONFIDENCE="call_min_confidence";
    static final double DEF_CALL_MIN_CONFIDENCE=30.0;
    static final double MIN_CALL_MIN_CONFIDENCE=0.0;
    static final double MAX_CALL_MIN_CONFIDENCE=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_call_min_confidence = new SettingsModelDoubleBounded(CFGKEY_CALL_MIN_CONFIDENCE, DEF_CALL_MIN_CONFIDENCE, MIN_CALL_MIN_CONFIDENCE, MAX_CALL_MIN_CONFIDENCE);
	// confidence threshold for emitting
	static final String CFGKEY_EMIT_MIN_CONFIDENCE="emit_min_confidence";
	static final double DEF_EMIT_MIN_CONFIDENCE=30.0;
    static final double MIN_EMIT_MIN_CONFIDENCE=0.0;
    static final double MAX_EMIT_MIN_CONFIDENCE=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_emit_min_confidence = new SettingsModelDoubleBounded(CFGKEY_EMIT_MIN_CONFIDENCE, DEF_EMIT_MIN_CONFIDENCE, MIN_EMIT_MIN_CONFIDENCE, MAX_EMIT_MIN_CONFIDENCE);
    // pcr error rate
	static final String CFGKEY_PCR_ERR = "pcr_error";
	static final double DEF_PCR_ERR=1e-4;
	static final double MIN_PCR_ERR = 0;
	static final double MAX_PCR_ERR=1;
	private final SettingsModelDoubleBounded m_pcr_error = new SettingsModelDoubleBounded(CFGKEY_PCR_ERR, DEF_PCR_ERR, MIN_PCR_ERR, MAX_PCR_ERR);
	// fraction of contamination -> reads to filter
	static final String CFGKEY_CONTAMINATION = "contamination";
	static final double DEF_CONTAMINATION=0;
	static final double MIN_CONTAMINATION=0;
	static final double MAX_CONTAMINATION=1;
	private final SettingsModelDoubleBounded m_contamination = new SettingsModelDoubleBounded(CFGKEY_CONTAMINATION, DEF_CONTAMINATION, MIN_CONTAMINATION, MAX_CONTAMINATION);
	// heterozygosity value -> how much 2 individuals differ
	static final String CFGKEY_HET="heterozygosity";
	static final double DEF_HET=0.001;
	static final double MIN_HET=0;
	static final double MAX_HET=1;
	private final SettingsModelDoubleBounded m_heterozygosity = new SettingsModelDoubleBounded(CFGKEY_HET, DEF_HET, MIN_HET, MAX_HET);
	// fraction of deletions -> threshold for calling
	static final String CFGKEY_MAX_DELETION_FRAC="max_deletion_fraction";
	static final double DEF_MAX_DELETION_FRAC=0.05;
	static final double MIN_MAX_DELETION_FRAC=0;
	static final double MAX_MAX_DELETION_FRAC=1.1;
	private final SettingsModelDoubleBounded m_max_deletion_fraction = new SettingsModelDoubleBounded(CFGKEY_MAX_DELETION_FRAC, DEF_MAX_DELETION_FRAC, MIN_MAX_DELETION_FRAC, MAX_MAX_DELETION_FRAC);
	// min base quality threshold
	static final String CFGKEY_MIN_BASE_QUAL="min_base_qual";
	static final int DEF_MIN_BASE_QUAL=17;
	static final int MIN_MIN_BASE_QUAL=0;
	static final int MAX_MIN_BASE_QUAL=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_base_qual = new SettingsModelIntegerBounded(CFGKEY_MIN_BASE_QUAL, DEF_MIN_BASE_QUAL, MIN_MIN_BASE_QUAL, MAX_MIN_BASE_QUAL);
	// indel heterozygosity
	static final String CFGKEY_INDEL_HET="indel_heterozygosity";
	static final double DEF_INDEL_HET=1.25e-4;
	static final double MIN_INDEL_HET=0;
	static final double MAX_INDEL_HET=1;
	private final SettingsModelDoubleBounded m_indel_heterozygosity = new SettingsModelDoubleBounded(CFGKEY_INDEL_HET, DEF_INDEL_HET, MIN_INDEL_HET, MAX_INDEL_HET);
	// minimum indel count for calling indel
	static final String CFGKEY_MIN_INDEL_CNT="min_indel_count";
	static final int DEF_MIN_INDEL_CNT=5;
	static final int MIN_MIN_INDEL_CNT=0;
	static final int MAX_MIN_INDEL_CNT=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_indel_count = new SettingsModelIntegerBounded(CFGKEY_MIN_INDEL_CNT, DEF_MIN_INDEL_CNT, MIN_MIN_INDEL_CNT, MAX_MIN_INDEL_CNT);
	// minimum indel fraction for calling indel -> only sites fullfilling this criteria will be considered in count threshold
	static final String CFGKEY_MIN_INDEL_FRAC="min_indel_frac";
	static final double DEF_MIN_INDEL_FRAC=0.25;
	static final double MIN_MIN_INDEL_FRAC=0;
	static final double MAX_MIN_INDEL_FRAC=1;
	private final SettingsModelDoubleBounded m_min_indel_frac = new SettingsModelDoubleBounded(CFGKEY_MIN_INDEL_FRAC, DEF_MIN_INDEL_FRAC, MIN_MIN_INDEL_FRAC, MAX_MIN_INDEL_FRAC);
	// indel gap open penalty
	static final String CFGKEY_GAP_OPEN_PEN = "gap_open_pen";
	static final int DEF_GAP_OPEN_PEN= 45;
	static final int MIN_GAP_OPEN_PEN=0;
	static final int MAX_GAP_OPEN_PEN=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_gap_open_pen = new SettingsModelIntegerBounded(CFGKEY_GAP_OPEN_PEN, DEF_GAP_OPEN_PEN, MIN_GAP_OPEN_PEN, MAX_GAP_OPEN_PEN);
	// indel continuation penalty
	static final String CFGKEY_GAP_CONT_PEN="gap_cont_pen";
	static final int DEF_GAP_CONT_PEN=10;
	static final int MIN_GAP_CONT_PEN=0;
	static final int MAX_GAP_CONT_PEN=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_gap_cont_pen = new SettingsModelIntegerBounded(CFGKEY_GAP_CONT_PEN, DEF_GAP_CONT_PEN, MIN_GAP_CONT_PEN, MAX_GAP_CONT_PEN);
	// malformed read filter MBQ = mismatching base and quals
	static final String CFGKEY_MBQ="mbq";
	static final boolean DEF_MBQ=true;
	private final SettingsModelBoolean m_mbq = new SettingsModelBoolean(CFGKEY_MBQ, DEF_MBQ);
	
	private int bam_index;
	private String outfile;
	private String dbsnp;

	//Network/Proxy options
//	public static final String CFGKEY_USEPROXY="useproxy";
//	public static final String CFGKEY_PROXYHOST="proxyhost";
//	public static final String CFGKEY_PROXYPORT="proxyport";
//	public static final String CFGKEY_USEPROXYAUTH="useproxyauth";
//	public static final String CFGKEY_PROXYUSER="proxyuser";
//	public static final String CFGKEY_PROXYPASSWORD="proxypassword";
	
//	private final SettingsModelBoolean m_useproxy = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXY, false);
//	private final SettingsModelString m_proxyhost = new SettingsModelString(
//			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYHOST,"");
//	private final SettingsModelString m_proxyport = new SettingsModelString(
//			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPORT,"");
//	private final SettingsModelBoolean m_useproxyauth = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXYAUTH, false);
//	private final SettingsModelString m_proxyuser = new SettingsModelString(
//			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYUSER,"");
//	private final SettingsModelString m_proxypassword = new SettingsModelString(
//			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPASSWORD,"");
	
    /**
     * Constructor for the node model.
     */
    protected GATKUnifiedGenotyperNodeModel() {
    
        // #in ports, #out ports
        super(OptionalPorts.createOPOs(1),OptionalPorts.createOPOs(1));

        addSetting(m_outfolder);
        addSetting(m_use_dbsnp);
        addSetting(m_dbsnp_file);
        addSetting(m_variant_type);
        addSetting(m_num_threads);
        
        addSetting(m_call_min_confidence);
        addSetting(m_emit_min_confidence);
        addSetting(m_pcr_error);
        addSetting(m_contamination);
        addSetting(m_heterozygosity);
        addSetting(m_max_deletion_fraction);
        addSetting(m_min_base_qual);
        addSetting(m_indel_heterozygosity);
        addSetting(m_min_indel_count);
        addSetting(m_min_indel_frac);
        addSetting(m_gap_open_pen);
        addSetting(m_gap_cont_pen);
        addSetting(m_mbq);
        
        addPrefPageSetting(m_dbsnp_file, IBISKNIMENodesPlugin.RES_DBSNP);
        
        //Proxy options
//    	m_proxyhost.setEnabled(false);
//    	m_proxyport.setEnabled(false);
//    	m_useproxyauth.setEnabled(false);
//    	m_proxyuser.setEnabled(false);
//    	m_proxypassword.setEnabled(false);
    }

    @Override
    protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {

    	//get and check input files
        ArrayList<String> inputfiles = new ArrayList<>();
        Iterator<DataRow> it = inData[0].iterator();
        while(it.hasNext()){
        	inputfiles.add(it.next().getCell(bam_index).toString());
        } 
        
        if(inputfiles.size()==0){
        	throw new InvalidSettingsException("No bam file available, something went wrong with the previous node!");
        }
        
        // check path(s) to bam file
        for(String in:inputfiles) {
        	if(!in.endsWith(".bam")) {
        		throw new InvalidSettingsException("Input file "+in+" is not in bam format!");
        	}
        	if(!Files.exists(Paths.get(in))){
            	throw new InvalidSettingsException("Path to input bam file: "+in+" does not exist");
            }
    		String index_file = IO.replaceFileExtension(in, "bai");

    		// check bam file index
    		if (Files.notExists(Paths.get(index_file)) && Files.notExists(Paths.get(in+".bai"))) {
    			throw new InvalidSettingsException("Missing BAM file index for "+in);
    		}
        }
        
        // create file names of vcf output  
        String type="";
        if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[0])){
        	type = "snps";
        }else if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[1])){
        	type = "indels";
        }else if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[2])){
        	type = "all_variants";
        }
        
        outfile = IO.replaceFileExtension(inputfiles.get(0), type + ".vcf");
        if(!CompatibilityChecker.inputFileNotOk(m_outfolder.getStringValue())) {
        	outfile = m_outfolder.getStringValue() + System.getProperty("file.separator")+ new File(inputfiles.get(0)).getName();
        	outfile = IO.replaceFileExtension(outfile, type + ".vcf");
        }
        
		//Enable proxy if needed
//		String proxyOptions = "";
//		if(m_useproxy.getBooleanValue()){
//			
//			proxyOptions += " -Dhttp.proxyHost=" + m_proxyhost.getStringValue();
//			proxyOptions += " -Dhttp.proxyPort=" + m_proxyport.getStringValue();
//			
//			if(m_useproxyauth.getBooleanValue()){
//				
//    			proxyOptions += " -Dhttp.proxyUser=" + m_proxyuser.getStringValue();
//    			proxyOptions += " -Dhttp.proxyPassword=" + m_proxypassword.getStringValue();
//			}
//			
//			proxyOptions += " ";
//		}
        
		String cmd="";
		cmd +=" -nt "+m_num_threads.getIntValue();
		
		for(String in:inputfiles){
			cmd+=" -I "+in;
		}
		
		cmd+=" -glm "+m_variant_type.getStringValue();
		
		if(m_use_dbsnp.getBooleanValue()){
        	cmd+= " -D " + dbsnp;
        }
		
		cmd+=" -stand_call_conf "+m_call_min_confidence.getDoubleValue();
		cmd+=" -stand_emit_conf "+m_emit_min_confidence.getDoubleValue();
		cmd+=" -pcr_error "+m_pcr_error.getDoubleValue();
		cmd+=" -contamination "+m_contamination.getDoubleValue();
		cmd+=" -hets "+m_heterozygosity.getDoubleValue();
		cmd+=" -deletions "+m_max_deletion_fraction.getDoubleValue();
		cmd+=" -mbq "+m_min_base_qual.getIntValue();
		cmd+=" -indelHeterozygosity "+m_indel_heterozygosity.getDoubleValue();
		cmd+=" -minIndelCnt "+ m_min_indel_count.getIntValue(); 
		cmd+=" -minIndelFrac "+m_min_indel_frac.getDoubleValue();
		cmd+=" -indelGOP "+ m_gap_open_pen.getIntValue();
		cmd+=" -indelGCP "+ m_gap_cont_pen.getIntValue();
		
		if(m_mbq.getBooleanValue()){
			cmd+=" --filter_mismatching_base_and_quals ";
		}
		
		return cmd;
    }

    @Override
    protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
    	bam_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("BAMCell")) {
    			bam_index = i;
    		}
    	}
		return (bam_index>-1);	
    }

	@Override
	protected String getCommandWalker() {
		return "UnifiedGenotyper";
	}

	@Override
	protected String getOutfile() {
		return this.outfile;
	}

	@Override
	protected DataType getOutColType() {
		return VCFCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		super.updatePrefs();
		dbsnp = m_dbsnp_file.getStringValue();
		
        if(m_use_dbsnp.getBooleanValue()) { 
        	if(CompatibilityChecker.inputFileNotOk(dbsnp)){
        		throw new InvalidSettingsException("Path to dbSNP file not given or incorrect!");
        	}
        	if(!Files.exists(Paths.get(dbsnp+".idx"))){
            	throw new InvalidSettingsException("dbSNP index file: "+dbsnp+".idx does not exist");
            }
        }
	}
}