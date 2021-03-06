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
package de.helmholtz_muenchen.ibis.ngs.snpsift;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;



/**
 * This is the model implementation of SnpSift.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class SnpSiftNodeModel extends HTExecutorNodeModel {
    
	static final String CFGKEY_SNPSIFT_BIN = "snpsift_bin";
	static final String CFGKEY_METHOD = "method";
	static final String DEF_METHOD = "";
	static final String DEF_FIELDS = "field1,field2";
	
	public enum SnpSiftTool {
		TSTV, FILTER, ANNOTATE, INTERVALS, DBNSFP, OTHER
	}
	static final LinkedHashMap<String, SnpSiftTool> NAME2TOOL = new LinkedHashMap<>();
	static {
		NAME2TOOL.put("Filter", SnpSiftTool.FILTER);
		NAME2TOOL.put("Annotate", SnpSiftTool.ANNOTATE);
		NAME2TOOL.put("TsTv", SnpSiftTool.TSTV);
		NAME2TOOL.put("Intervals", SnpSiftTool.INTERVALS);
		NAME2TOOL.put("dbNSFP", SnpSiftTool.DBNSFP);
		NAME2TOOL.put("Other", SnpSiftTool.OTHER);
	}
	
	//Filter
	public static final String CFGKEY_FILTERSTRING = "filterstring";
	
	//Annotate
	public static final String CFGKEY_ANNID = "annid";
	public static final String CFGKEY_ANNINFO = "anninfo";
	public static final String CFGKEY_ANNVCFDB = "annvcfdb";
	public static final String CFGKEY_ANN_OPT = "ann_opt_field";
	
	//Intervals
	public static final String CFGKEY_INTERX = "interx";
	public static final String CFGKEY_INTERBED = "interbed";
	
	//dbNFSP
	public static final String CFGKEY_DBNSFP = "dbnsfp";
	public static final String CFGKEY_DBNSFPFFIELDS = "dbnsfpFields";
	public static final String CFGKEY_DBNSFPFFIELDSALL = "dbnsfpFieldsall";
	public static final String CFGKEY_DBNSFP_OPT = "dbnsfp_opt_field";
	
	//Other
	public static final String CFGKEY_OTHER_CMD = "other_cmd";
	public static final String CFGKEY_OTHER_OUT = "other_out";
	
	/**Setting Models**/
	private final SettingsModelString m_snpsift_bin = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_SNPSIFT_BIN,"");
	private final SettingsModelString m_method = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_METHOD,SnpSiftNodeModel.DEF_METHOD);
	
	/**Filter**/
	private final SettingsModelString m_filterstring = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_FILTERSTRING,"");

	/**Annotate**/
	private final SettingsModelOptionalString m_anninfo = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_ANNINFO,DEF_FIELDS,false);
	private final SettingsModelBoolean m_annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
	private final SettingsModelString m_annvcfdb = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
	private final SettingsModelOptionalString m_ann_opt = new SettingsModelOptionalString(SnpSiftNodeModel.CFGKEY_ANN_OPT, "", false);
	
	/**Intervals**/
	private final SettingsModelString m_interbed = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INTERBED,"");
	private final SettingsModelBoolean m_interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

	/**dbNSFP**/
	private final SettingsModelString m_dbnsfp = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
	private final SettingsModelOptionalString m_dbnsfpfields = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,DEF_FIELDS,false);
	private final SettingsModelBoolean m_dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);
	private final SettingsModelOptionalString m_dbnsfp_opt = new SettingsModelOptionalString(SnpSiftNodeModel.CFGKEY_DBNSFP_OPT, "", false);

	/*Other*/
	private final SettingsModelString m_other_cmd = new SettingsModelString(SnpSiftNodeModel.CFGKEY_OTHER_CMD,"");
	private final SettingsModelString m_other_out = new SettingsModelString(SnpSiftNodeModel.CFGKEY_OTHER_OUT,"");
	
	private int vcf_index;
	private String snpsift_bin;
	private DataType outType = FileCell.TYPE;

	public static final String OUT_COL = "snpSift_result";
	
    protected SnpSiftNodeModel() {
        super(1, 1, 1);
        
		addSetting(m_method);
		addSetting(m_snpsift_bin);
		addSetting(m_filterstring);
		addSetting(m_annid);
		addSetting(m_anninfo);
		addSetting(m_annvcfdb);
		addSetting(m_ann_opt);
		addSetting(m_interbed);
		addSetting(m_interx);
		addSetting(m_dbnsfp);
		addSetting(m_dbnsfpfields);
		addSetting(m_dbnsfpfieldsall);
		addSetting(m_dbnsfp_opt);
		addSetting(m_other_cmd);
		addSetting(m_other_out);
		
		addPrefPageSetting(m_snpsift_bin, IBISKNIMENodesPlugin.SNPSIFT);
		
		m_other_cmd.setEnabled(false);
		m_other_out.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
        String stdOutFile = "";
        ArrayList<String> command = new ArrayList<String>();
        command.add("java");
    	command.add("-jar");
    	command.add(snpsift_bin);
        
        String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();

        SnpSiftTool tool = NAME2TOOL.get(m_method.getStringValue());
        
        switch(tool) {
        case TSTV:
        	command.add("tstv");
        	command.add(vcf_infile);
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_tstv.txt");
        	break;
        case FILTER:
        	command.add("filter");
        	command.add("--file");
        	command.add(vcf_infile);
        	command.add(m_filterstring.getStringValue());
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_filtered.vcf");
        	break;
        case ANNOTATE:
        	command.add("annotate");
        	if(m_annid.getBooleanValue()){
        		command.add("-noInfo");
        	}
        	if(m_anninfo.isActive()){
        		command.add("-info");
        		command.add(m_anninfo.getStringValue());
        	}
        	if(m_ann_opt.getStringValue().length() > 1) {
        		String [] opts = m_ann_opt.getStringValue().split(" ");
        		for(String s:opts) {
        			command.add(s.trim());
        		}
        	}
        	command.add(IO.processFilePath(m_annvcfdb.getStringValue()));
        	command.add(vcf_infile);
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_annotated.vcf");
        	break;
        case INTERVALS:
        	command.add("intervals");
        	command.add("-i");
        	command.add(vcf_infile);
        	if(m_interx.getBooleanValue()){
        		command.add("-x");
        	}
        	command.add(IO.processFilePath(m_interbed.getStringValue()));
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_intervals.vcf");
        	break;
        case DBNSFP:
        	command.add("dbnsfp");
        	if(m_dbnsfpfieldsall.getBooleanValue()){
        		command.add("-a");
        	}
        	if(m_dbnsfpfields.isActive()){
        		command.add("-f");
        		command.add(m_dbnsfpfields.getStringValue());
        	}
        	if(m_dbnsfp_opt.getStringValue().length() > 1) {
        		String [] opts = m_dbnsfp_opt.getStringValue().split(" ");
        		for(String s:opts) {
        			command.add(s.trim());
        		}
        	}
        	command.add("-db");
        	command.add(IO.processFilePath(m_dbnsfp.getStringValue()));
        	command.add(vcf_infile);
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_dbnfsp.vcf");
        	break;
        case OTHER:
        	String [] cmd = m_other_cmd.getStringValue().split(" ");
        	for(String c:cmd) {
        		command.add(c.trim());
        	}
        	command.add(vcf_infile);
        	stdOutFile = IO.processFilePath(m_other_out.getStringValue());
        	break;
		default:
			break;
        }

        
    	/**Execute**/
        String [] cmd = new String [command.size()];
    	for(int i=0; i < cmd.length; i++) {
    		cmd[i] = command.get(i);
    	}
    	String lockFile = stdOutFile + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(cmd, stdOutFile, exec, new File(lockFile),stdOutFile);
	
    	//Create Output Table
    	
    	if(stdOutFile.endsWith(".vcf")) {
    		outType = VCFCell.TYPE;
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL, outType).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(stdOutFile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	super.updatePrefs();
    	snpsift_bin = IO.processFilePath(m_snpsift_bin.getStringValue());
    	
    	if(SnpSiftNodeDialog.getUseMainInputColBool()){
    		vcf_index = inSpecs[0].findColumnIndex(SnpSiftNodeDialog.getMainInputCol1());
    		if(!inSpecs[0].getColumnSpec(vcf_index).getType().toString().equals("VCFCell")){
    			vcf_index = -1;
    		}
    	} else {
    		vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
    	}
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}
    	
    	outType = FileCell.TYPE;
    	SnpSiftTool tool = 	NAME2TOOL.get(m_method.getStringValue());
    	switch(tool) {
    	case FILTER:
    		if(m_filterstring.getStringValue().equals("")) {
    			throw new InvalidSettingsException("No filter string defined!");
    		}
    		outType = VCFCell.TYPE;
    		break;
    	case ANNOTATE:
    		if(CompatibilityChecker.inputFileNotOk(IO.processFilePath(m_annvcfdb.getStringValue()))) {
    			throw new InvalidSettingsException("Given annotation database invalid!");
    		}
    		outType = VCFCell.TYPE;
    		break;
    	case INTERVALS:
    		if(CompatibilityChecker.inputFileNotOk(IO.processFilePath(m_interbed.getStringValue()))) {
    			throw new InvalidSettingsException("Interval file invalid!");
    		}
    		outType = VCFCell.TYPE;
    		break;
    	case DBNSFP:
    		if(CompatibilityChecker.inputFileNotOk(IO.processFilePath(m_dbnsfp.getStringValue()))) {
    			throw new InvalidSettingsException("dbNSFP file invalid!");
    		}
    		if(CompatibilityChecker.inputFileNotOk(IO.processFilePath(m_dbnsfp.getStringValue() + ".tbi"))) {
    			throw new InvalidSettingsException("dbNSFP file invalid!");
    		}
    		outType = VCFCell.TYPE;
    		break;
    	case OTHER:
    		if(m_other_cmd.getStringValue().equals("") || m_other_out.getStringValue().equals("")) {
    			throw new InvalidSettingsException("Both, command and output file have to be defined!");
    		}
    		if(Files.notExists(Paths.get(new File(IO.processFilePath(m_other_out.getStringValue())).getParent()))) {
    			throw new InvalidSettingsException("Directory for output file does not exist!");
    		}
    		if(IO.processFilePath(m_other_out.getStringValue()).endsWith(".vcf")) {
    			outType = VCFCell.TYPE;
    		}
    		break;
    	default:
    		break;
    	}
    	
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL, outType).createSpec()})};
    }
}

