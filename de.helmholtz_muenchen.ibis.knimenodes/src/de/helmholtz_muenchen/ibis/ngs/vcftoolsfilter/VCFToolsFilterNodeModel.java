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
package de.helmholtz_muenchen.ibis.ngs.vcftoolsfilter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
/**
 * This is the model implementation of LOFFilter.
 * 
 *
 * @author Tim Jeske
 */
public class VCFToolsFilterNodeModel extends HTExecutorNodeModel {
    
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFToolsFilterNodeModel.class);
	
    
    static final String CFGKEY_VCFTOOLS = "vcf_tools";
    private final SettingsModelString m_vcf_tools = new SettingsModelString(CFGKEY_VCFTOOLS,"-");
    
    static final String CFGKEY_FILL_AN_AC = "fill_an_ac";
    final SettingsModelBoolean m_fill_an_ac = new SettingsModelBoolean(CFGKEY_FILL_AN_AC,false);
    
    
    // output options
    //static final String CFGKEY_OUTFOLDER = "outfolder";
    //final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
    
    //genotype filter
    static final String CFGKEY_FILTER_BY_DP = "filter_by_DP";
    final SettingsModelBoolean m_filter_by_DP = new SettingsModelBoolean(CFGKEY_FILTER_BY_DP,false);
    
    static final String CFGKEY_DP_THRESHOLD = "DP_threshold";
    final SettingsModelInteger m_DP_threshold = new SettingsModelInteger(CFGKEY_DP_THRESHOLD,8);
    
    static final String CFGKEY_FILTER_BY_AD = "filter_by_AD";
    final SettingsModelBoolean m_filter_by_AD = new SettingsModelBoolean(CFGKEY_FILTER_BY_AD,false);
    
    static final String CFGKEY_AD_THRESHOLD = "AD_threshold";
    final SettingsModelInteger m_AD_threshold = new SettingsModelInteger(CFGKEY_AD_THRESHOLD,8);
    
    static final String CFGKEY_FILTER_BY_GQ = "filter_by_GQ";
    final SettingsModelBoolean m_filter_by_GQ = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ,false);
    
    static final String CFGKEY_GQ_THRESHOLD = "GQ_threshold";
    final SettingsModelIntegerBounded m_GQ_threshold = new SettingsModelIntegerBounded(CFGKEY_GQ_THRESHOLD,20,0,99);
    
    //variant filter
    static final String CFGKEY_FILTER_PASS = "filter_pass";
    final SettingsModelBoolean m_filter_pass = new SettingsModelBoolean(CFGKEY_FILTER_PASS,false);
    
    static final String CFGKEY_FILTER_CALL_RATE = "filter_callRate";
    final SettingsModelBoolean m_filter_by_callRate = new SettingsModelBoolean(CFGKEY_FILTER_CALL_RATE,false);
    
    static final String CFGKEY_CALL_RATE_THRESHOLD = "callRate_threshold";
    final SettingsModelDoubleBounded m_callRate_threshold = new SettingsModelDoubleBounded(CFGKEY_CALL_RATE_THRESHOLD, 0.88,0.0,1.0);
    
    static final String CFGKEY_FILTER_BY_GQ_MEAN = "filter_by_GQ_MEAN";
    final SettingsModelBoolean m_filter_by_GQ_MEAN = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ_MEAN,false);
    
    static final String CFGKEY_GQ_MEAN_THRESHOLD = "GQ_mean_threshold";
    final SettingsModelIntegerBounded m_GQ_mean_threshold = new SettingsModelIntegerBounded(CFGKEY_GQ_MEAN_THRESHOLD,35,0,99);
    
    
	//output col names
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
	private int vcf_index;
	private String vcftools_bin;
	
    protected VCFToolsFilterNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    	
    	//addSetting(m_outfolder);
    	addSetting(m_filter_by_DP);
    	addSetting(m_filter_by_AD);
    	addSetting(m_DP_threshold);
    	addSetting(m_filter_by_GQ);
    	addSetting(m_GQ_threshold);
    	addSetting(m_filter_pass);
    	addSetting(m_filter_by_callRate);
    	addSetting(m_callRate_threshold);
    	addSetting(m_filter_by_GQ_MEAN);
    	addSetting(m_GQ_mean_threshold);
    	addSetting(m_vcf_tools);
        addSetting(m_fill_an_ac);
        
        addPrefPageSetting(m_vcf_tools, IBISKNIMENodesPlugin.VCFTOOLS);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	//check input file
    	String infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String outfile = infile;
    	if(m_filter_by_DP.getBooleanValue() || m_filter_by_GQ.getBooleanValue() || m_filter_by_AD.getBooleanValue()) {
    		outfile = filterGenotypes(infile, exec);
    		infile = outfile;
    	}
    	if(m_filter_pass.getBooleanValue() || m_filter_by_GQ_MEAN.getBooleanValue() || m_filter_by_callRate.getBooleanValue()) {
    		outfile = filterVariants(infile, exec);
    		infile = outfile;
    	}
    	if(m_fill_an_ac.getBooleanValue()) {
    		outfile = postprocess(infile,exec);
    	}
    	
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	return new BufferedDataTable[]{outTable};
    }

    private String filterGenotypes(String infile, ExecutionContext exec) throws Exception {
    	
    	boolean is_gz = infile.endsWith(".gz");
    	String outfile;
    	if(is_gz) {
    		outfile = infile.replace(".vcf.gz", ".genotype_filtered");
    	} else {
    		outfile = IO.replaceFileExtension(infile, ".genotype_filtered");
    	}
    	
    	if(m_filter_by_GQ.getBooleanValue() || m_filter_by_DP.getBooleanValue()) {
    		ArrayList<String> cmd = new ArrayList<>();
    		cmd.add(vcftools_bin);
    		if(is_gz) {
    			cmd.add("--gzvcf");
    		} else {
    			cmd.add("--vcf");
    		}
    		cmd.add(infile);
    		cmd.add("--out");
    		cmd.add(outfile);
    		cmd.add("--minDP");
    		cmd.add(m_DP_threshold.getIntValue()+"");
    		cmd.add("--minGQ");
    		cmd.add(m_GQ_threshold.getIntValue()+"");
    		cmd.add("--recode");
    		cmd.add("--recode-INFO-all");
    		
    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
    		outfile = outfile+".recode.vcf";
    		File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);
    		super.executeCommand(cmd_array, outfile, exec, lockFile);
    	}
    	
    	if(m_filter_by_AD.getBooleanValue()) {
    		
    		if(infile.endsWith(".gz")) {
    			LOGGER.warn("AD filter can only process unzipped vcf files!");
    			return outfile;
    		}
    		
    		outfile = infile.replace(".vcf",".AD_filtered.vcf");
    		
    		ArrayList<String> cmd = new ArrayList<>();
			cmd.add("perl");
			cmd.add(IO.getScriptPath()+"scripts/perl/filterByAD.pl");
			cmd.add(infile);
			cmd.add(m_AD_threshold.getIntValue()+"");
			cmd.add("./.");
			
			String[] cmd_array = new String[cmd.size()];
			for (int i = 0; i < cmd.size(); i++) {
				cmd_array[i] = cmd.get(i);
			}

			File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);

			super.executeCommand(cmd_array, outfile, exec, lockFile);
    		
    	}
    	
    	return outfile;
    }
    
    private String filterVariants(String infile, ExecutionContext exec) throws Exception {
    	String outfile = "";
    	
    	if(m_filter_pass.getBooleanValue() || m_filter_by_callRate.getBooleanValue()) {
    		boolean is_gz = infile.endsWith(".gz");
    		if(is_gz) {
        		outfile = infile.replace(".vcf.gz", ".variant_filtered");
        	} else {
        		outfile = IO.replaceFileExtension(infile, ".variant_filtered");
        	}
    		
    		ArrayList<String> cmd = new ArrayList<>();
    		cmd.add(vcftools_bin);
    		if(is_gz) {
    			cmd.add("--gzvcf");
    		} else {
    			cmd.add("--vcf");
    		}
    		cmd.add(infile);
    		cmd.add("--out");
    		cmd.add(outfile);
    		if(m_filter_pass.getBooleanValue()) {
    			cmd.add("--remove-filtered-all");
    		}
    		if(m_filter_by_callRate.getBooleanValue()) {
    			cmd.add("--max-missing");
    			cmd.add(m_callRate_threshold.getDoubleValue()+"");
    		}
    		
    		cmd.add("--recode");
    		cmd.add("--recode-INFO-all");
    		
    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
    		outfile = outfile+".recode.vcf";
    		File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);
    		
    		super.executeCommand(cmd_array, outfile, exec, lockFile);
    		infile = outfile;
    	}
    	
    	if(m_filter_by_GQ_MEAN.getBooleanValue()) {
    		
    		if(infile.endsWith(".gz")) {
    			LOGGER.warn("Mean GQ filter can only process unzipped vcf files!");
    			return outfile;
    		}
    		
    		outfile = infile.replace(".vcf",".GQFiltered.vcf");
    		
    		ArrayList<String> cmd = new ArrayList<>();
			cmd.add("perl");
			cmd.add(IO.getScriptPath()+"scripts/perl/filterGQMean.pl");
			cmd.add(infile);
			cmd.add(m_GQ_mean_threshold.getIntValue()+"");
			cmd.add("0");
			
			String[] cmd_array = new String[cmd.size()];
			for (int i = 0; i < cmd.size(); i++) {
				cmd_array[i] = cmd.get(i);
			}

			File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);

			super.executeCommand(cmd_array, outfile, exec, lockFile);
    	}
    	return outfile;
    }
    
    private String postprocess(String infile, ExecutionContext exec) throws Exception {
    	boolean is_gz = infile.endsWith(".gz");
    	String outfile;
    	if(is_gz) {
    		outfile = infile.replace(".vcf.gz", ".adjusted_AN_AC.vcf");
    	} else {
    		outfile = infile.replace(".vcf", ".adjusted_AN_AC.vcf");
    	}
    	
    	String cmd = "";
    	if(is_gz) {
			cmd = "zcat";
		} else {
			cmd = "cat";
		}
		
		cmd += " " + infile;
		
    	if(m_fill_an_ac.getBooleanValue()) {
    		
    		String vcftools = vcftools_bin;
    		int pos = vcftools.lastIndexOf(System.getProperty("file.separator"));
    		vcftools = vcftools.substring(0,pos)+ System.getProperty("file.separator")+"fill-an-ac";

    		cmd += " | " + vcftools;
    	}

    	File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);	
    	cmd += " > " + outfile;
    	super.executeCommand(new String[]{"/usr/bin/bash","-c",cmd}, outfile ,exec,lockFile);
    	
    	return outfile;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	super.updatePrefs();
    	vcftools_bin = IO.processFilePath(m_vcf_tools.getStringValue());
    	
    	boolean a = m_filter_by_DP.getBooleanValue();
    	boolean b = m_filter_by_GQ.getBooleanValue();
    	boolean c = m_filter_by_callRate.getBooleanValue();
    	boolean d = m_fill_an_ac.getBooleanValue();
    	boolean f = m_filter_pass.getBooleanValue();
    	vcf_index = -1;
    	
    	if(a || b || c || d || f) {
    		String vcftools_warning = CheckUtils.checkSourceFile(vcftools_bin);
        	if(vcftools_warning != null) {
        		setWarningMessage(vcftools_warning);
        	}
    	}
    	
    	for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
    
    	
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }
}

