package de.helmholtz_muenchen.ibis.ngs.bowtie2;

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

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * This is the model implementation of Bowtie2.
 * 
 *
 * @author 
 */
public class Bowtie2NodeModel extends HTExecutorNodeModel {
    
    	public static final String CFGKEY_INSTALLPATH = "installpath";
    	public static final String CFGKEY_REFSEQFILE = "refseqfile";
    	public static final String CFGKEY_NOAUTO = "noauto";
    	public static final String CFGKEY_PACKED = "packed";
    	public static final String CFGKEY_BMAX = "bmax";
    	public static final String CFGKEY_DCV = "dcv";
    	public static final String CFGKEY_NODC = "nodc";
    	public static final String CFGKEY_OFFRATE = "offrate";
    	public static final String CFGKEY_FTABCHARS = "ftabchars";
    	public static final String CFGKEY_USECUTOFF = "usecutoff";

    	public static final String CFGKEY_USESKIP = "useskip";
    	public static final String CFGKEY_SKIP = "skip";
    	public static final String CFGKEY_USEUPTO = "useupto";
    	public static final String CFGKEY_UPTO = "upto";
    	public static final String CFGKEY_TRIM5 = "trim5";
    	public static final String CFGKEY_TRIM3 = "trim3";
    	public static final String CFGKEY_QUALS = "quals";
    	public static final String CFGKEY_USEPRESET = "usepreset";
    	public static final String CFGKEY_PRESET = "preset";
    	public static final String CFGKEY_N = "n";
    	public static final String CFGKEY_L = "l";
    	public static final String CFGKEY_I1 = "i1";
    	public static final String CFGKEY_I2 = "i2";
    	public static final String CFGKEY_NCEIL1 = "nceil1";
    	public static final String CFGKEY_NCEIL2 = "nceil2";
    	public static final String CFGKEY_DPAD = "dpad";
    	public static final String CFGKEY_GBAR = "gbar";
    	public static final String CFGKEY_IGNOREQUALS = "ignorequals";
    	public static final String CFGKEY_NOFW = "nofw";
    	public static final String CFGKEY_NORC = "norc";
    	public static final String CFGKEY_ALIGNMENTTYPE = "alignmenttype";
    	public static final String CFGKEY_MA = "ma";
    	public static final String CFGKEY_MP = "mp";
    	public static final String CFGKEY_NP = "np";
    	public static final String CFGKEY_RDG1 = "rdg1";
    	public static final String CFGKEY_RDG2 = "rdg2";
    	public static final String CFGKEY_RFG1 = "rfg1";
    	public static final String CFGKEY_RFG2 = "rfg2";
    	public static final String CFGKEY_SCOREMIN1 = "scoremin1";
    	public static final String CFGKEY_SCOREMIN2 = "scoremin2";
    	public static final String CFGKEY_REPORTING1 = "reporting1";
    	public static final String CFGKEY_REPORTING2 = "reporting2";
    	public static final String CFGKEY_D = "d";
    	public static final String CFGKEY_R = "r";
    	public static final String CFGKEY_MININS = "minins";
    	public static final String CFGKEY_MAXINS = "maxins";
    	public static final String CFGKEY_FF = "ff";
    	public static final String CFGKEY_NOMIXED = "nomixed";
    	public static final String CFGKEY_NODISCORDANT = "nodiscordant";
    	public static final String CFGKEY_NODOVETAIL = "nodovetail";
    	public static final String CFGKEY_NOCONTAIN = "nocontain";
    	public static final String CFGKEY_NOOVERLAP = "nooverlap";
    	public static final String CFGKEY_THREADS = "threads";
    	public static final String CFGKEY_RECORDER = "recorder";
    	public static final String CFGKEY_MM = "mm";
    	public static final String CFGKEY_QCFILTER = "qcfilter";


    	public static final int DEFAULT_BMAX = 4;
    	public static final int DEFAULT_DCV = 1024;
    	public static final int DEFAULT_OFFRATE = 5;
    	public static final int DEFAULT_FTABCHARS = 10;
    	public static final int DEFAULT_CUTOFF = Integer.MAX_VALUE;
    	public static final int DEFAULT_SKIP = 0;
    	public static final int DEFAULT_UPTO = Integer.MAX_VALUE;
    	public static final int DEFAULT_TRIM5 = 0;
    	public static final int DEFAULT_TRIM3 = 0;
    	public static final int DEFAULT_N = 0;
    	public static final int DEFAULT_L = 22;
    	public static final double DEFAULT_I1 = 1;
    	public static final double DEFAULT_I2 = 1.15;
    	public static final double DEFAULT_NCEIL1 = 0;
    	public static final double DEFAULT_NCEIL2 = 0.15;
    	public static final int DEFAULT_DPAD = 15;
    	public static final int DEFAULT_GBAR = 4;
    	public static final int DEFAULT_MA = 2;
    	public static final int DEFAULT_MP = 6;
    	public static final int DEFAULT_NP = 1;
    	public static final int DEFAULT_RDG1 = 5;
    	public static final int DEFAULT_RDG2 = 3;
    	public static final int DEFAULT_RFG1 = 5;
    	public static final int DEFAULT_RFG2 = 3;
    	public static final double DEFAULT_SCOREMIN1 = -0.6;
    	public static final double DEFAULT_SCOREMIN2 = -0.6;
    	public static final int DEFAULT_REPORTING2 = 1;
    	public static final int DEFAULT_D = 15;
    	public static final int DEFAULT_R = 2;
    	public static final int DEFAULT_MININS = 0;
    	public static final int DEFAULT_MAXINS = 500;
    	public static final int DEFAULT_THREADS = 4;

    	private final SettingsModelString m_installpath = new SettingsModelString(CFGKEY_INSTALLPATH,"");
    	private final SettingsModelString m_refseqfile = new SettingsModelString(CFGKEY_REFSEQFILE,"");
    	private final SettingsModelBoolean m_noauto = new SettingsModelBoolean(CFGKEY_NOAUTO, true);
    	private final SettingsModelBoolean m_packed = new SettingsModelBoolean(CFGKEY_PACKED, false);
    	private final SettingsModelIntegerBounded m_bmax = new SettingsModelIntegerBounded(CFGKEY_BMAX,DEFAULT_BMAX,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_dcv = new SettingsModelIntegerBounded(CFGKEY_DCV,DEFAULT_DCV,1,4096);
    	private final SettingsModelBoolean m_nodc = new SettingsModelBoolean(CFGKEY_NODC, false);
    	private final SettingsModelIntegerBounded m_offrate = new SettingsModelIntegerBounded(CFGKEY_OFFRATE,DEFAULT_OFFRATE,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_ftabchars = new SettingsModelIntegerBounded(CFGKEY_FTABCHARS,DEFAULT_FTABCHARS,1,Integer.MAX_VALUE);
    	private final SettingsModelBoolean m_useskip = new SettingsModelBoolean(CFGKEY_USESKIP, false);
    	private final SettingsModelIntegerBounded m_skip = new SettingsModelIntegerBounded(CFGKEY_SKIP,DEFAULT_SKIP,0,Integer.MAX_VALUE);
    	private final SettingsModelBoolean m_useupto = new SettingsModelBoolean(CFGKEY_USEUPTO, false);
    	private final SettingsModelIntegerBounded m_upto = new SettingsModelIntegerBounded(CFGKEY_UPTO,DEFAULT_UPTO,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_trim5 = new SettingsModelIntegerBounded(CFGKEY_TRIM5,DEFAULT_TRIM5,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_trim3 = new SettingsModelIntegerBounded(CFGKEY_TRIM3,DEFAULT_TRIM3,0,Integer.MAX_VALUE);
    	private final SettingsModelString m_quals = new SettingsModelString(CFGKEY_QUALS,"");
    	private final SettingsModelBoolean m_usepreset = new SettingsModelBoolean(CFGKEY_USEPRESET, true);
    	private final SettingsModelString m_preset = new SettingsModelString(CFGKEY_PRESET,"");
    	private final SettingsModelIntegerBounded m_n = new SettingsModelIntegerBounded(CFGKEY_N,DEFAULT_N,0,1);
    	private final SettingsModelIntegerBounded m_l = new SettingsModelIntegerBounded(CFGKEY_L,DEFAULT_L,4,31);
    	private final SettingsModelDoubleBounded m_i1 = new SettingsModelDoubleBounded(CFGKEY_I1,DEFAULT_I1,0,Double.MAX_VALUE);
    	private final SettingsModelDoubleBounded m_i2 = new SettingsModelDoubleBounded(CFGKEY_I2,DEFAULT_I2,0,Double.MAX_VALUE);
    	private final SettingsModelDoubleBounded m_nceil1 = new SettingsModelDoubleBounded(CFGKEY_NCEIL1,DEFAULT_NCEIL1,0,Double.MAX_VALUE);
    	private final SettingsModelDoubleBounded m_nceil2 = new SettingsModelDoubleBounded(CFGKEY_NCEIL2,DEFAULT_NCEIL2,0,Double.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_dpad = new SettingsModelIntegerBounded(CFGKEY_DPAD,DEFAULT_DPAD,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_gbar = new SettingsModelIntegerBounded(CFGKEY_GBAR,DEFAULT_GBAR,0,Integer.MAX_VALUE);
    	private final SettingsModelBoolean m_ignorequals = new SettingsModelBoolean(CFGKEY_IGNOREQUALS, false);
    	private final SettingsModelBoolean m_nofw = new SettingsModelBoolean(CFGKEY_NOFW, false);
    	private final SettingsModelBoolean m_norc = new SettingsModelBoolean(CFGKEY_NORC, false);
    	private final SettingsModelString m_alignmenttype = new SettingsModelString(CFGKEY_ALIGNMENTTYPE,"");
    	private final SettingsModelIntegerBounded m_ma = new SettingsModelIntegerBounded(CFGKEY_MA,DEFAULT_MA,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_mp = new SettingsModelIntegerBounded(CFGKEY_MP,DEFAULT_MP,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_np = new SettingsModelIntegerBounded(CFGKEY_NP,DEFAULT_NP,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_rdg1 = new SettingsModelIntegerBounded(CFGKEY_RDG1,DEFAULT_RDG1,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_rdg2 = new SettingsModelIntegerBounded(CFGKEY_RDG2,DEFAULT_RDG2,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_rfg1 = new SettingsModelIntegerBounded(CFGKEY_RFG1,DEFAULT_RFG1,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_rfg2 = new SettingsModelIntegerBounded(CFGKEY_RFG2,DEFAULT_RFG2,0,Integer.MAX_VALUE);
    	private final SettingsModelDoubleBounded m_scoremin1 = new SettingsModelDoubleBounded(CFGKEY_SCOREMIN1,DEFAULT_SCOREMIN1,-Double.MAX_VALUE,Double.MAX_VALUE);
    	private final SettingsModelDoubleBounded m_scoremin2 = new SettingsModelDoubleBounded(CFGKEY_SCOREMIN2,DEFAULT_SCOREMIN2,-Double.MAX_VALUE,Double.MAX_VALUE);
    	private final SettingsModelString m_reporting1 = new SettingsModelString(CFGKEY_REPORTING1,"");
    	private final SettingsModelIntegerBounded m_reporting2 = new SettingsModelIntegerBounded(CFGKEY_REPORTING2,DEFAULT_REPORTING2,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_d = new SettingsModelIntegerBounded(CFGKEY_D,DEFAULT_D,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_r = new SettingsModelIntegerBounded(CFGKEY_R,DEFAULT_R,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_minins = new SettingsModelIntegerBounded(CFGKEY_MININS,DEFAULT_MININS,0,Integer.MAX_VALUE);
    	private final SettingsModelIntegerBounded m_maxins = new SettingsModelIntegerBounded(CFGKEY_MAXINS,DEFAULT_MAXINS,0,Integer.MAX_VALUE);
    	private final SettingsModelString m_ff = new SettingsModelString(CFGKEY_FF,"");
    	private final SettingsModelBoolean m_nomixed = new SettingsModelBoolean(CFGKEY_NOMIXED, false);
    	private final SettingsModelBoolean m_nodiscordant = new SettingsModelBoolean(CFGKEY_NODISCORDANT, false);
    	private final SettingsModelBoolean m_nodovetail = new SettingsModelBoolean(CFGKEY_NODOVETAIL, false);
    	private final SettingsModelBoolean m_nocontain = new SettingsModelBoolean(CFGKEY_NOCONTAIN, false);
    	private final SettingsModelBoolean m_nooverlap = new SettingsModelBoolean(CFGKEY_NOOVERLAP, false);
    	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(CFGKEY_THREADS,DEFAULT_THREADS,1,Integer.MAX_VALUE);
    	private final SettingsModelBoolean m_recorder = new SettingsModelBoolean(CFGKEY_RECORDER, false);
    	private final SettingsModelBoolean m_mm = new SettingsModelBoolean(CFGKEY_MM, true);
    	private final SettingsModelBoolean m_qcfilter = new SettingsModelBoolean(CFGKEY_QCFILTER, false);

		
    	//The Output Col Names
    	public static final String OUT_COL1 = "Path2SAMFile";
    	
    	private static String readType = "";
    	
    /**
     * Constructor for the node model.
     */
    protected Bowtie2NodeModel() {
    	
        super(1, 1);
        
        m_packed.setEnabled(false);
        m_bmax.setEnabled(false);
        m_dcv.setEnabled(false);
        m_skip.setEnabled(false);
    	m_upto.setEnabled(false);
    	m_n.setEnabled(false);
    	m_l.setEnabled(false);
    	m_i1.setEnabled(false);
    	m_i2.setEnabled(false);
    	m_preset.setStringValue("sensitive");
    	m_alignmenttype.setStringValue("entire read must align (no clipping)");
    	m_ma.setEnabled(false);
    	m_reporting2.setEnabled(false);
    	m_d.setEnabled(false);
    	m_r.setEnabled(false);
    	m_minins.setEnabled(false);
    	m_maxins.setEnabled(false);
    	m_ff.setEnabled(false);
    	m_nomixed.setEnabled(false);
    	m_nodiscordant.setEnabled(false);
    	m_nodovetail.setEnabled(false);
    	m_nocontain.setEnabled(false);
    	m_nooverlap.setEnabled(false);
    	
    	
    	addSetting(m_d);
    	addSetting(m_dcv);
    	addSetting(m_dpad);
    	addSetting(m_ff);
    	addSetting(m_ftabchars);
    	addSetting(m_gbar);
    	addSetting(m_i1);
    	addSetting(m_i2);
    	addSetting(m_ignorequals);
    	addSetting(m_installpath);
    	addSetting(m_l);
    	addSetting(m_ma);
    	addSetting(m_maxins);
    	addSetting(m_minins);
    	addSetting(m_mm);
    	addSetting(m_mp);
    	addSetting(m_n);
    	addSetting(m_nceil1);
    	addSetting(m_nceil2);
    	addSetting(m_noauto);
    	addSetting(m_nocontain);
    	addSetting(m_nodc);
    	addSetting(m_nodiscordant);
    	addSetting(m_nodovetail);
    	addSetting(m_nofw);
    	addSetting(m_nomixed);
    	addSetting(m_nooverlap);
    	addSetting(m_norc);
    	addSetting(m_np);
    	addSetting(m_offrate);
    	addSetting(m_packed);
    	addSetting(m_preset);
    	addSetting(m_qcfilter);
    	addSetting(m_quals);
    	addSetting(m_r);
    	addSetting(m_rdg1);
    	addSetting(m_rdg2);
    	addSetting(m_recorder);
    	addSetting(m_refseqfile);
    	addSetting(m_reporting1);
    	addSetting(m_reporting2);
    	addSetting(m_rfg1);
    	addSetting(m_rfg2);
    	addSetting(m_scoremin1);
    	addSetting(m_scoremin2);
    	addSetting(m_skip);
    	addSetting(m_threads);
    	addSetting(m_trim3);
    	addSetting(m_trim5);
    	addSetting(m_upto);
    	addSetting(m_usepreset);
    	addSetting(m_useskip);
    	addSetting(m_useupto);

    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//Prepare File
    	String path2refFile 	 = m_refseqfile.getStringValue();
    	String path2readFile 	 = inData[0].iterator().next().getCell(0).toString();
    	String path2readFile2 	 = "NA";
    	if(readType.equals("paired-end")){
    		path2readFile2 = inData[0].iterator().next().getCell(1).toString();
    	}

    	//Prepare Outfile Name
    	String basePath 	= path2readFile.substring(0,path2readFile.lastIndexOf("."));
    	String path2outfile = basePath+".sam";
    	
    	/**
    	 * INDEX REFGENOME
    	 */

    	createIndex(exec);
   
    	//execute
    	generateAlignment(exec,path2refFile,path2readFile,path2readFile2,path2outfile);
     	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, SAMCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2outfile)};
    	
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
    	

    	CompatibilityChecker CC = new CompatibilityChecker();
    	readType = CC.getReadType(inSpecs, 0);
    	if(CC.getWarningStatus()){
    		setWarningMessage(CC.getWarningMessages());
    	}
    	if(readType.equals("paired-end")){
    		m_minins.setEnabled(true);
        	m_maxins.setEnabled(true);
        	m_ff.setEnabled(true);
        	m_nomixed.setEnabled(true);
        	m_nodiscordant.setEnabled(true);
        	m_nodovetail.setEnabled(true);
        	m_nocontain.setEnabled(true);
        	m_nooverlap.setEnabled(true);
    	}else{
    		m_minins.setEnabled(false);
        	m_maxins.setEnabled(false);
        	m_ff.setEnabled(false);
        	m_nomixed.setEnabled(false);
        	m_nodiscordant.setEnabled(false);
        	m_nodovetail.setEnabled(false);
        	m_nocontain.setEnabled(false);
        	m_nooverlap.setEnabled(false);
    	}

    	
    	if(m_refseqfile.getStringValue().length() > 1) {
    		if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	            throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
	    	}
    	}

    	String path2baseFileName = "";
    	if(m_refseqfile.getStringValue().isEmpty()){
    		throw new InvalidSettingsException("Reference Genome is missing!");
    	}else{
    		path2baseFileName = m_refseqfile.getStringValue().substring(0,m_refseqfile.getStringValue().lastIndexOf("."));
    	}
    	
    	
    	//Check for bowtie2-align
    	Boolean f1 = !new File(path2baseFileName + ".1.bt2").exists();
    	Boolean f2 = !new File(path2baseFileName + ".2.bt2").exists();
    	Boolean f3 = !new File(path2baseFileName + ".3.bt2").exists();
    	Boolean f4 = !new File(path2baseFileName + ".4.bt2").exists();
    	Boolean f5 = !new File(path2baseFileName + ".rev.1.bt2").exists();
    	Boolean f6 = !new File(path2baseFileName + ".rev.2.bt2").exists();
    	if(f1 && f2 && f3 && f4 && f5 && f6) {
        	String bowtie2_align = m_installpath.getStringValue().replaceAll("bowtie2$", "bowtie2-build");
        	if(!new File(bowtie2_align).exists()){
        		throw new InvalidSettingsException("bowtie2-build binary is missing. However it should be located at "+bowtie2_align);
        	}
    	}
    	
    	

    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, SAMCell.TYPE).createSpec()})};
    }

    /**
     * Build bowtie2 index
     * @param exec
     * @throws Exception
     */
    protected void createIndex(ExecutionContext exec) throws Exception{
    	
    	// Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
    	
    	ArrayList<String> command = new ArrayList<String>();
    	String path2refFile = m_refseqfile.getStringValue();
    	String path2baseFileName = path2refFile.substring(0,path2refFile.lastIndexOf("."));
    	String bowtieBuild = m_installpath.getStringValue().replace("bowtie2", "bowtie2-build");

    	command.add(bowtieBuild);
   
    	if(!m_noauto.getBooleanValue()) {
    		command.add("--noauto");
    		if(m_packed.getBooleanValue()) {
    			command.add("--packed");
    		}
    		command.add("--bmaxdivn " + m_bmax.getIntValue());
    		command.add("--dcv " + m_dcv.getIntValue());
    		}
    	
    	if(m_nodc.getBooleanValue()) {
    		command.add("--nodc");
    	}
    	command.add("--offrate " + m_offrate.getIntValue());
    	command.add("--ftabchars " + m_ftabchars.getIntValue());

    	command.add(path2refFile);
    	command.add(path2baseFileName);

    	/**Execute**/
    	String lockFile = path2refFile + "_bowtie2_index" + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));        	
    }
    
    /**
     * Generates bowtie2 alignment
     * @param exec
     * @param path2refFile
     * @param path2readFile
     * @param path2readFile2
     * @param path2outfile
     * @throws Exception
     */
    protected void generateAlignment(ExecutionContext exec, String path2refFile, String path2readFile, String path2readFile2, String path2outfile) throws Exception{
    	ArrayList<String> command = new ArrayList<String>();
    	String bowtieAlign = m_installpath.getStringValue();
    	String path2baseFileName = path2refFile.substring(0,path2refFile.lastIndexOf("."));


    	command.add(bowtieAlign);
    	
    	if(m_useskip.getBooleanValue()) {
    		command.add("--skip " + m_skip.getIntValue());
    	}
    	if(m_useupto.getBooleanValue()) {
    		command.add("--upto " + m_upto.getIntValue());
    	}
    	if(m_trim5.getIntValue() > 0) {
    		command.add("--trim5 " + m_trim5.getIntValue()); 
    	}
    	if(m_trim3.getIntValue() > 0) {
    		command.add("--trim3 " + m_trim3.getIntValue()); 
    	}
    	if(m_quals.getStringValue().equals("Phred+64")) {
    		command.add("--phred64");
    	} else if(m_quals.getStringValue().equals("Encoded as space-delimited integers")) {
    		command.add("--int-quals");
    	}
    	
    	if(m_usepreset.getBooleanValue()) {
    		String ps = m_preset.getStringValue();
    		String preset = "";
    		if(ps.equals("very-fast")) {
    			preset = " --very-fast";
    		} else if(ps.equals("fast")) {
    			preset = " --fast";
    		} else if(ps.equals("sensitive")) {
    			preset = " --sensitive";
    		} else if(ps.equals("very-sensitive")) {
    			preset = " --very-sensitive";
    		}
    		if(m_alignmenttype.getStringValue().equals("local alignment (ends might be soft clipped)")) {
    			preset += "-local";
    		}
    		command.add(preset);
    		
    	} else {
    		command.add("-N " + m_n.getIntValue());
    		command.add("-L " + m_l.getIntValue());
    		command.add("-i S," + m_i1.getDoubleValue()+","+m_i2.getDoubleValue());
    		if(m_alignmenttype.getStringValue().equals("local alignment (ends might be soft clipped)")) {
    			command.add("--local");
    		}
    	}
    	command.add("--n-ceil " + m_nceil1.getDoubleValue() + "," + m_nceil2.getDoubleValue());
    	command.add("--dpad " + m_dpad.getIntValue());
    	command.add("--gbar " + m_gbar.getIntValue());
    	if(m_ignorequals.getBooleanValue()) {
    		command.add(" --ignore-quals");
    	}
    	if(m_nofw.getBooleanValue()) {
    		command.add("--nofw");
    	}
    	if(m_norc.getBooleanValue()) {
    		command.add("--norc");
    	}
    	if(m_ma.isEnabled()) {
    		command.add("--ma " + m_ma.getIntValue());
    	}
    	command.add("--mp " + m_mp.getIntValue());
    	command.add("--np " + m_np.getIntValue());
    	command.add("--rdg " + m_rdg1.getIntValue() + "," + m_rdg2.getIntValue());
    	command.add("--rfg " + m_rfg1.getIntValue() + "," + m_rfg2.getIntValue());
    	String scoremin = "";
    	scoremin = "--score-min L," + m_scoremin1.getDoubleValue() + "," + m_scoremin2.getDoubleValue();
    	if(m_alignmenttype.getStringValue().equals("local alignment (ends might be soft clipped)")) {
    		scoremin = "--score-min G," + m_scoremin1.getDoubleValue() + "," + m_scoremin2.getDoubleValue();
    	}
    	command.add(scoremin);
    	String reporting = "";
    	if(m_reporting1.getStringValue().equals("Report up to <int> alns per read; MAPQ not meaningful:")) {
    		reporting = "-k " + m_reporting2.getIntValue();
    	} else if(m_reporting1.getStringValue().equals("Report all alignments; very slow, MAPQ not meaningful")) {
    		reporting = "--all";
    	}
    	command.add(reporting);
    	
    	if(m_d.isEnabled()) {
    		command.add("-D " + m_d.getIntValue());
    	}
    	if(m_r.isEnabled()) {
    		command.add("-R " + m_r.getIntValue());
    	}
    	
    	if(readType.equals("paired-end")) {
    		command.add("--minins " + m_minins.getIntValue());
    		command.add("--maxins " + m_maxins.getIntValue());
	    	if(m_ff.getStringValue().equals("reverse/forward")) {
	    		command.add("--rf");
	    	} else if(m_ff.getStringValue().equals("forward/forward")) {
	    		command.add("--ff");
	    	}
	    	
	    	if(m_nomixed.getBooleanValue()) {
	    		command.add("--no-mixed");
	    	}
	    	if(m_nodiscordant.getBooleanValue()) {
	    		command.add("--no-discordant");
	    	}
	    	if(m_nodovetail.getBooleanValue()) {
	    		command.add("--no-dovetail");
	    	}
	    	if(m_nocontain.getBooleanValue()) {
	    		command.add("--no-contain");
	    	}
	    	if(m_nooverlap.getBooleanValue()) {
	    		command.add("--no-overlap");
	    	}
	    	
    	}
    	command.add("--threads " + m_threads.getIntValue());
    	if(m_recorder.getBooleanValue()) {
    		command.add("--reorder");
    	}
    	if(m_mm.getBooleanValue()) {
    		command.add("--mm");
    	}
    	if(m_qcfilter.getBooleanValue()) {
    		command.add("--qc-filter");
    	}
    	  
    	String reads = "";
    	if(readType.equals("single-end")) {
    		reads = "-U " + path2readFile;
    	} else {
    		reads = "-1 " + path2readFile + " -2 " + path2readFile2;
    	}
    	
    	command.add("-x " + path2baseFileName);
    	command.add(reads);
    	command.add("-S " + path2outfile);

    	/**Execute**/
    	String lockFile = path2outfile + "_bowtie2" + SuccessfulRunChecker.LOCK_ENDING;
    	String stdOut 	= path2outfile+".stdOut";
    	String stdErr 	= path2outfile+".stdErr";
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile),stdOut,stdErr);  
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

