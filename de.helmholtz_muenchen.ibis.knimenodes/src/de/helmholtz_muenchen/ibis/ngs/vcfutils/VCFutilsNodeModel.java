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


package de.helmholtz_muenchen.ibis.ngs.vcfutils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.knime.core.data.DataTableSpec;
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
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of VCFutils.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class VCFutilsNodeModel extends HTExecutorNodeModel {
	
	public static final String CFGKEY_VCFFILE = "vcffile";
	public static final String CFGKEY_UTILITY = "utility";
	public static final String CFGKEY_VCF = "vcf";
	public static final String CFGKEY_REFVCFFILE = "refvcffile";
	public static final String CFGKEY_QSTATSUSEREF = "qstatusref";
	public static final String CFGKEY_SNPFILE = "snpfile";
	public static final String CFGKEY_HAPMAPFILE = "hapmapfile";
	public static final String CFGKEY_MINRMS = "minrms";
	public static final String CFGKEY_MINREADDEPTH = "minreaddepth";
	public static final String CFGKEY_MAXREADDEPTH = "maxreaddepth";
	public static final String CFGKEY_MINALTBASE = "minaltbases";
	public static final String CFGKEY_GAPFILTER = "gapfilter";
	public static final String CFGKEY_ADJACENTGAPS = "adjacentgaps";
	public static final String CFGKEY_STRANDPVAL = "strandpval";
	public static final String CFGKEY_BASEQPVAL = "basqpval";
	public static final String CFGKEY_MAPQPVAL = "mapqpval";
	public static final String CFGKEY_ENDDISTPVAL = "enddistpval";
	public static final String CFGKEY_HWEPVAL = "hwepval";
	public static final String CFGKEY_PRINTFILTERED = "printfiltered";
	public static final String CFGKEY_INDELFW = "indelfw";

	public static final int DEFAULT_MINRMS = 10;
	public static final int DEFAULT_MINREADDEPTH = 2;
	public static final int DEFAULT_MAXREADDEPTH = 10000000;
	public static final int DEFAULT_MINALTBASE = 2;
	public static final int DEFAULT_GAPFILTER = 3;
	public static final int DEFAULT_ADJACENTGAPS = 10;
	public static final double DEFAULT_STRANDPVAL = 0.0001;
	public static final int DEFAULT_BASEQPVAL = 100;
	public static final double DEFAULT_MAPQPVAL = 0;
	public static final double DEFAULT_ENDDISTPVAL = 0.0001;
	public static final double DEFAULT_HWEPVAL = 0.0001;
	public static final int DEFAULT_INDELFW = 5;

	private final SettingsModelString m_vcffile = new SettingsModelString(VCFutilsNodeModel.CFGKEY_VCFFILE,"");
	private final SettingsModelString m_utility = new SettingsModelString(VCFutilsNodeModel.CFGKEY_UTILITY,"subsam");
	private final SettingsModelString m_vcf = new SettingsModelString(VCFutilsNodeModel.CFGKEY_VCF,"");
	private final SettingsModelString m_refvcffile = new SettingsModelString(VCFutilsNodeModel.CFGKEY_REFVCFFILE,"");
	private final SettingsModelBoolean m_qstatusref = new SettingsModelBoolean(VCFutilsNodeModel.CFGKEY_QSTATSUSEREF, false);
	private final SettingsModelString m_snpfile = new SettingsModelString(VCFutilsNodeModel.CFGKEY_SNPFILE,"");
	private final SettingsModelString m_hapmapfile = new SettingsModelString(VCFutilsNodeModel.CFGKEY_HAPMAPFILE,"");
	private final SettingsModelIntegerBounded m_minrms = new SettingsModelIntegerBounded(CFGKEY_MINRMS,DEFAULT_MINRMS,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minreaddepth = new SettingsModelIntegerBounded(CFGKEY_MINREADDEPTH,DEFAULT_MINREADDEPTH,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxreaddepth = new SettingsModelIntegerBounded(CFGKEY_MAXREADDEPTH,DEFAULT_MAXREADDEPTH,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minaltbase = new SettingsModelIntegerBounded(CFGKEY_MINALTBASE,DEFAULT_MINALTBASE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_gapfilter = new SettingsModelIntegerBounded(CFGKEY_GAPFILTER,DEFAULT_GAPFILTER,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_adjacentgaps = new SettingsModelIntegerBounded(CFGKEY_ADJACENTGAPS,DEFAULT_ADJACENTGAPS,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_strandpval = new SettingsModelDoubleBounded(CFGKEY_STRANDPVAL,DEFAULT_STRANDPVAL,0,Double.MAX_VALUE);
	private final SettingsModelIntegerBounded m_baseqpval = new SettingsModelIntegerBounded(CFGKEY_BASEQPVAL,DEFAULT_BASEQPVAL,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_mapqpval = new SettingsModelDoubleBounded(CFGKEY_MAPQPVAL,DEFAULT_MAPQPVAL,0,Double.MAX_VALUE);
	private final SettingsModelDoubleBounded m_enddistpval = new SettingsModelDoubleBounded(CFGKEY_ENDDISTPVAL,DEFAULT_ENDDISTPVAL,0,Double.MAX_VALUE);
	private final SettingsModelDoubleBounded m_hwepval = new SettingsModelDoubleBounded(CFGKEY_HWEPVAL,DEFAULT_HWEPVAL,0,Double.MAX_VALUE);	
	private final SettingsModelBoolean m_printfiltered = new SettingsModelBoolean(CFGKEY_PRINTFILTERED, false);
	private final SettingsModelIntegerBounded m_indelfw = new SettingsModelIntegerBounded(CFGKEY_INDELFW,DEFAULT_INDELFW,0,Integer.MAX_VALUE);

	public static boolean optionalPort=false;	
	
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFutilsNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected VCFutilsNodeModel() {
    	
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(0), 2);
        
        readType = "paired-end";
        
    	addSetting(m_adjacentgaps);
    	addSetting(m_baseqpval);
    	addSetting(m_enddistpval);
    	addSetting(m_gapfilter);
    	addSetting(m_hapmapfile);
    	addSetting(m_hwepval);
    	addSetting(m_indelfw);
    	addSetting(m_mapqpval);
    	addSetting(m_maxreaddepth);
    	addSetting(m_minaltbase);
    	addSetting(m_minreaddepth);
    	addSetting(m_minrms);
    	addSetting(m_printfiltered);
    	addSetting(m_qstatusref);
    	addSetting(m_refvcffile);
    	addSetting(m_snpfile);
    	addSetting(m_strandpval);
    	addSetting(m_utility);
    	addSetting(m_vcf);
    	addSetting(m_vcffile);
        
        m_qstatusref.setEnabled(false);
        m_refvcffile.setEnabled(false);
        m_snpfile.setEnabled(false);
        m_hapmapfile.setEnabled(false);
        m_minrms.setEnabled(false);
        m_minreaddepth.setEnabled(false);
        m_maxreaddepth.setEnabled(false);
        m_minaltbase.setEnabled(false);
        m_gapfilter.setEnabled(false);
        m_adjacentgaps.setEnabled(false);
        m_strandpval.setEnabled(false);
    	m_baseqpval.setEnabled(false);
    	m_mapqpval.setEnabled(false);
    	m_enddistpval.setEnabled(false);
    	m_hwepval.setEnabled(false);
    	m_printfiltered.setEnabled(false);
    	m_indelfw.setEnabled(false);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String fle = "";
    	if(optionalPort){
    		int colIndex2 = 1;
        	if(VCFutilsNodeDialog.getUseMainInputColBool()){
        		colIndex2 = inData[0].getDataTableSpec().findColumnIndex(VCFutilsNodeDialog.getMainInputCol2());
        	}
    		fle = inData[0].iterator().next().getCell(colIndex2).toString();
    	}
    	else{
    		if(m_vcffile.isEnabled()) {
    			fle = m_vcffile.getStringValue();
    		}	
    		else{
    			fle = m_snpfile.getStringValue();
    		}
    	}
    	
    	String path2vcf = "";
    	if(optionalPort){
        	int colIndex1 = 0;
        	if(VCFutilsNodeDialog.getUseMainInputColBool()){
        		colIndex1 = inData[0].getDataTableSpec().findColumnIndex(VCFutilsNodeDialog.getMainInputCol1());
        	}
    		String path2bcftools = inData[0].iterator().next().getCell(colIndex1).toString();
    		path2vcf = path2bcftools+"vcfutils.pl";
    	}
    	else{
    		path2vcf = m_vcf.getStringValue();
    	}
    	String utility = m_utility.getStringValue();
    	String path2vcffile = "";
    	String minrms = "-Q "+m_minrms.getIntValue();
    	String minreaddepth = "-d " + m_minreaddepth.getIntValue();
    	String maxreaddepth = "-D " + m_maxreaddepth.getIntValue();
    	String minaltbases = "-a " + m_minaltbase.getIntValue();
    	String gapfilter = "-w " + m_gapfilter.getIntValue();
    	String adjacentgaps = "-W " + m_adjacentgaps.getIntValue();
    	String strandpval = "-1 " + m_strandpval.getDoubleValue();
    	String baseqpval = "-2 1.0E-" + m_baseqpval.getIntValue();
    	String mapqpval = "-3 " + m_mapqpval.getDoubleValue();
    	String enddistpval = "-4 " + m_enddistpval.getDoubleValue();
    	String hwepval = "-e " + m_hwepval.getDoubleValue();
    	String printfiltered = "";
    	String indelfilteringwindow = "-l " + m_indelfw.getIntValue();
    	
    	if(m_printfiltered.getBooleanValue()) {
    		printfiltered="-p ";
    	}
    	
    	if(m_vcffile.isEnabled()) {
    		path2vcffile = m_vcffile.getStringValue();
    	}else if(optionalPort){
    		path2vcffile=fle;
    	}
    	
    	
    	ArrayList<String> command = new ArrayList<String>();
    	String outfile = "";
    	command.add(path2vcf);
    	
    // "subsam", "listsam", "fillac", "qstats", "hapmap2vcf", "ucscsnp2vcf", "varFilter", "vcf2fq"
    	if(utility.equals("subsam")) {
    // vcfutils.pl subsam <in.vcf> [samples]
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+"_subsam.vcf";
    		command.add(utility);
    		command.add(path2vcffile);
    	} else if(utility.equals("listsam")) {
    // vcfutils.pl listsam <in.vcf>
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+"_listsam.txt";
    		command.add(utility);
    		command.add(path2vcffile);
    	} else if(utility.equals("fillac")) {
    // vcfutils.pl fillac <in.vcf>
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+"_fillac.vcf";
    		command.add(utility);
    		command.add(path2vcffile);
    	} else if(utility.equals("qstats")) {
    // vcfutils.pl qstats [-r ref.vcf] <in.vcf>
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+"_qstats.txt";
    		command.add(utility);
    		if(m_qstatusref.getBooleanValue()) {
    			command.add("-r " + m_refvcffile.getStringValue());
    		}
    		command.add(path2vcffile);
    	} else if(utility.equals("hapmap2vcf")) {
    // vcfutils.pl <in.ucsc.snp> <in.hapmap>
    		String p2hapmapfile = m_hapmapfile.getStringValue();
    		outfile = p2hapmapfile.substring(0, p2hapmapfile.lastIndexOf("."))+".vcf";
    		command.add(m_snpfile.getStringValue());
    		command.add(p2hapmapfile);
    	} else if(utility.equals("ucscsnp2vcf")) {
    // vcfutils.pl <in.ucsc.snp>
    		String p2snpfile = m_snpfile.getStringValue();
    		outfile = p2snpfile.substring(0, p2snpfile.lastIndexOf("."))+".vcf";
    		command.add(p2snpfile);
    	} else if(utility.equals("varFilter")) {
    // vcfutils.pl varFilter [options] <in.vcf>
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+"_filtered.vcf";
    		command.add(utility);
    		command.add(minrms);
    		command.add(minreaddepth);
    		command.add(maxreaddepth);
    		command.add(minaltbases);
    		command.add(gapfilter);
    		command.add(adjacentgaps);
    		command.add(strandpval);
    		command.add(baseqpval);
    		command.add(mapqpval);
    		command.add(enddistpval);
    		command.add(hwepval);
    		command.add(printfiltered);
    		command.add(path2vcffile);
    	} else if(utility.equals("vcf2fq")) {
    // vcfutils.pl vcf2fq [options] <all-site.vcf>
    		outfile = path2vcffile.substring(0, path2vcffile.lastIndexOf("."))+".fq";
    		command.add(utility);
    		command.add(minrms);
    		command.add(minreaddepth);
    		command.add(maxreaddepth);
    		command.add(indelfilteringwindow);
    		command.add(path2vcffile);
    	}
    	
    	/**
    	 * Execute
    	 */
    	String lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{command.toString()}, outfile, exec, new File(lockFile),outfile);
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,outfile);
  	
        return new BufferedDataTable[]{};
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
    	
    	String vcfinfile="";
    	//Check OptionalInputPort
		try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
			//Get INFILE Values
			if(getAvailableInputFlowVariables().containsKey("BCFTOOLSOUT")){
				//BCFTOOLS already executed
				 vcfinfile=getAvailableInputFlowVariables().get("BCFTOOLSOUT").getStringValue();
				String suffix=vcfinfile.substring(vcfinfile.lastIndexOf("."));
				if(!suffix.equals(".vcf")){
					throw new InvalidSettingsException("VCFUtils requires VCF input file."+suffix+" format found. For 'hapmap2vcf' and 'ucscsnp2vcf' delete the node connection.");
				}
			}
				//Format OK
		
			String[] colnames = inSpecs[0].getColumnNames();
			if(VCFutilsNodeDialog.getUseMainInputColBool()){
        		if(VCFutilsNodeDialog.getMainInputCol1().equals("Path2Bcftools") && VCFutilsNodeDialog.getMainInputCol1().equals("Path2BcftoolsOutput")){
        			m_vcffile.setEnabled(false);
    				m_vcf.setEnabled(false);
    				m_snpfile.setEnabled(false);
    				m_hapmapfile.setEnabled(false);
        		} else {
        			throw new InvalidSettingsException("Expected first selected column's name to be 'Path2Bcftools' and second 'Path2BcftoolsOutput'.");
        		}
        	}
			else if(colnames[0].equals("Path2Bcftools") && colnames[1].equals("Path2BcftoolsOutput")){
				m_vcffile.setEnabled(false);
				m_vcf.setEnabled(false);
				m_snpfile.setEnabled(false);
				m_hapmapfile.setEnabled(false);
			}
			else{
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			}
    	}catch(NullPointerException e){
    		vcfinfile=m_vcffile.getStringValue();
    		m_vcffile.setEnabled(true);
    		m_vcf.setEnabled(true);
    		if(m_utility.getStringValue().equals("hapmap2vcf")){
        		m_snpfile.setEnabled(true);
        		m_hapmapfile.setEnabled(true);
    		}
    		if(m_utility.getStringValue().equals("ucscsnp2vcf")){
    			m_snpfile.setEnabled(true);
    		}
    	}
    	
    	if(optionalPort && (m_utility.getStringValue().equals("hapmap2vcf") || m_utility.getStringValue().equals("ucscsnp2vcf") )){
    		throw new InvalidSettingsException("The utility you selected is not available when another node is connected to this node's in-port!");
    	}
    	
    	
    	
    	if(m_vcffile.isEnabled()) {
    		String filept = m_vcffile.getStringValue();
    		if(filept.length() > 1) {
	    		String fileTail = filept.substring(filept.length()-3,filept.length());
	    		if(!fileTail.equals("vcf") && !fileTail.equals("VCF")) {
	    			throw new InvalidSettingsException("The vcf file needs to have the file extension 'vcf'!");
	    		}
    		}
    	}
    	
    	if(m_refvcffile.isEnabled()) {
    		String filept = m_refvcffile.getStringValue();
    		if(filept.length() > 1) {
	    		String fileTail = filept.substring(filept.length()-3,filept.length());
	    		if(!fileTail.equals("vcf") && !fileTail.equals("VCF")) {
	    			throw new InvalidSettingsException("The reference vcf file needs to have the file extension 'vcf'!");
	    		}
    		}
    	}
    	
    	if(m_snpfile.isEnabled()) {
    		String filept = m_snpfile.getStringValue();
    		if(filept.length() > 1) {
	    		String fileTail = filept.substring(filept.length()-3,filept.length());
	    		if(!fileTail.equals("snp") && !fileTail.equals("SNP")) {
	    			throw new InvalidSettingsException("The snp file needs to have the file extension 'snp'!");
	    		}
    		}
    	}
    	
    	if(m_hapmapfile.isEnabled()) {
    		String filept = m_hapmapfile.getStringValue();
    		if(filept.length() > 1) {
	    		String fileTail = filept.substring(filept.length()-6,filept.length());
	    		if(!fileTail.equals("hapmap") && !fileTail.equals("HAPMAP")) {
	    			throw new InvalidSettingsException("The hapmap file needs to have the file extension 'hapmap'!");
	    		}
    		}
    	}
    	
		/**
		 * Check if Outfile exists
		 */		
		String path2outputfile2="";
		if(vcfinfile.length()>0){
		    // "subsam", "listsam", "fillac", "qstats", "hapmap2vcf", "ucscsnp2vcf", "varFilter", "vcf2fq"
	    	if(m_utility.getStringValue().equals("subsam")) {
	    // vcfutils.pl subsam <in.vcf> [samples]
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+"_subsam.vcf";
	    	} else if(m_utility.getStringValue().equals("listsam")) {
	    // vcfutils.pl listsam <in.vcf>
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+"_listsam.txt";
	    	} else if(m_utility.getStringValue().equals("fillac")) {
	    // vcfutils.pl fillac <in.vcf>
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+"_fillac.vcf";
	    	} else if(m_utility.getStringValue().equals("qstats")) {
	    // vcfutils.pl qstats [-r ref.vcf] <in.vcf>
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+"_qstats.txt";
	    	} else if(m_utility.getStringValue().equals("hapmap2vcf")) {
	    // vcfutils.pl <in.ucsc.snp> <in.hapmap>
	    		String p2hapmapfile = m_hapmapfile.getStringValue();
	    		 path2outputfile2 = p2hapmapfile.substring(0, p2hapmapfile.lastIndexOf("."))+".vcf";
	    	} else if(m_utility.getStringValue().equals("ucscsnp2vcf")) {
	    // vcfutils.pl <in.ucsc.snp>
	    		String p2snpfile = m_snpfile.getStringValue();
	    		 path2outputfile2 = p2snpfile.substring(0, p2snpfile.lastIndexOf("."))+".vcf";
	    	} else if(m_utility.getStringValue().equals("varFilter")) {
	    // vcfutils.pl varFilter [options] <in.vcf>
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+"_filtered.vcf";
	    	} else if(m_utility.getStringValue().equals("vcf2fq")) {
	    // vcfutils.pl vcf2fq [options] <all-site.vcf>
	    		 path2outputfile2 = vcfinfile.substring(0, vcfinfile.lastIndexOf("."))+".fq";	    		
		}						
		}
		System.out.println("############");
		System.out.println(m_utility.getStringValue());
    	String outfile=path2outputfile2;
    	System.out.println(path2outputfile2);
    	File outpath = new File(outfile);
    	if(outpath.exists()){
    		setWarningMessage("Outfile "+outfile + " already exists ! Please rename or move to other directory.");
    	}
    	
    	
        return new DataTableSpec[]{};
    }

//	/**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void saveSettingsTo(final NodeSettingsWO settings) {
//    	
//    	/** added for HTE **/
//    	super.saveSettingsTo(settings);
//    	
//    	m_adjacentgaps.saveSettingsTo(settings);
//    	m_baseqpval.saveSettingsTo(settings);
//    	m_enddistpval.saveSettingsTo(settings);
//    	m_gapfilter.saveSettingsTo(settings);
//    	m_hapmapfile.saveSettingsTo(settings);
//    	m_hwepval.saveSettingsTo(settings);
//    	m_indelfw.saveSettingsTo(settings);
//    	m_mapqpval.saveSettingsTo(settings);
//    	m_maxreaddepth.saveSettingsTo(settings);
//    	m_minaltbase.saveSettingsTo(settings);
//    	m_minreaddepth.saveSettingsTo(settings);
//    	m_minrms.saveSettingsTo(settings);
//    	m_printfiltered.saveSettingsTo(settings);
//    	m_qstatusref.saveSettingsTo(settings);
//    	m_refvcffile.saveSettingsTo(settings);
//    	m_snpfile.saveSettingsTo(settings);
//    	m_strandpval.saveSettingsTo(settings);
//    	m_utility.saveSettingsTo(settings);
//    	m_vcf.saveSettingsTo(settings);
//    	m_vcffile.saveSettingsTo(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//    	
//    	/** added for HTE **/
//    	super.loadValidatedSettingsFrom(settings);
//    	
//    	m_adjacentgaps.loadSettingsFrom(settings);
//    	m_baseqpval.loadSettingsFrom(settings);
//    	m_enddistpval.loadSettingsFrom(settings);
//    	m_gapfilter.loadSettingsFrom(settings);
//    	m_hapmapfile.loadSettingsFrom(settings);
//    	m_hwepval.loadSettingsFrom(settings);
//    	m_indelfw.loadSettingsFrom(settings);
//    	m_mapqpval.loadSettingsFrom(settings);
//    	m_maxreaddepth.loadSettingsFrom(settings);
//    	m_minaltbase.loadSettingsFrom(settings);
//    	m_minreaddepth.loadSettingsFrom(settings);
//    	m_minrms.loadSettingsFrom(settings);
//    	m_printfiltered.loadSettingsFrom(settings);
//    	m_qstatusref.loadSettingsFrom(settings);
//    	m_refvcffile.loadSettingsFrom(settings);
//    	m_snpfile.loadSettingsFrom(settings);
//    	m_strandpval.loadSettingsFrom(settings);
//    	m_utility.loadSettingsFrom(settings);
//    	m_vcf.loadSettingsFrom(settings);
//    	m_vcffile.loadSettingsFrom(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void validateSettings(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//    	
//    	/** added for HTE **/
//    	super.validateSettings(settings);
//    	
//    	m_adjacentgaps.validateSettings(settings);
//    	m_baseqpval.validateSettings(settings);
//    	m_enddistpval.validateSettings(settings);
//    	m_gapfilter.validateSettings(settings);
//    	m_hapmapfile.validateSettings(settings);
//    	m_hwepval.validateSettings(settings);
//    	m_indelfw.validateSettings(settings);
//    	m_mapqpval.validateSettings(settings);
//    	m_maxreaddepth.validateSettings(settings);
//    	m_minaltbase.validateSettings(settings);
//    	m_minreaddepth.validateSettings(settings);
//    	m_minrms.validateSettings(settings);
//    	m_printfiltered.validateSettings(settings);
//    	m_qstatusref.validateSettings(settings);
//    	m_refvcffile.validateSettings(settings);
//    	m_snpfile.validateSettings(settings);
//    	m_strandpval.validateSettings(settings);
//    	m_utility.validateSettings(settings);
//    	m_vcf.validateSettings(settings);
//    	m_vcffile.validateSettings(settings);
//    }
    
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

