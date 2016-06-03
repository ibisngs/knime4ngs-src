package de.helmholtz_muenchen.ibis.ngs.snpsift;

import java.io.File;
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
	static final String DEF_METHOD = "Filter";
	
	public enum SnpSiftTool {
		TSTV, FILTER, ANNOTATE, INTERVALS, ANNOTATE_DBSNP
	}
	static final LinkedHashMap<String, SnpSiftTool> NAME2TOOL = new LinkedHashMap<>();
	static {
		NAME2TOOL.put("Filter", SnpSiftTool.FILTER);
		NAME2TOOL.put("Annotate", SnpSiftTool.ANNOTATE);
		NAME2TOOL.put("TsTv", SnpSiftTool.TSTV);
		NAME2TOOL.put("Intervals", SnpSiftTool.INTERVALS);
		NAME2TOOL.put("Annotate with dbsnfp", SnpSiftTool.ANNOTATE_DBSNP);
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
	
	//dbnsfp
	public static final String CFGKEY_DBNSFP = "dbnsfp";
	public static final String CFGKEY_DBNSFPFFIELDS = "dbnsfpFields";
	public static final String CFGKEY_DBNSFPFFIELDSALL = "dbnsfpFieldsall";
	
	
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
			SnpSiftNodeModel.CFGKEY_ANNINFO,"",false);
	private final SettingsModelBoolean m_annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
	private final SettingsModelString m_annvcfdb = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
	private final SettingsModelOptionalString m_ann_opt = new SettingsModelOptionalString(SnpSiftNodeModel.CFGKEY_ANN_OPT, "", false);
	
	/**Intervals**/
	private final SettingsModelString m_interbed = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INTERBED,"");
	private final SettingsModelBoolean m_interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

	/**dbnsfp**/
	private final SettingsModelString m_dbnsfp = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
	private final SettingsModelOptionalString m_dbnsfpfields = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,"",false);
	private final SettingsModelBoolean m_dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);

	private int vcf_index;
	private DataType outType = FileCell.TYPE;

	public static final String OUT_COL = "snpSift_result";
	
    protected SnpSiftNodeModel() {
        super(1, 1);
        
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
        String stdOutFile = "";
        ArrayList<String> command = new ArrayList<String>();
        command.add("java");
    	command.add("-jar");
    	command.add(m_snpsift_bin.getStringValue());
        
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
            	command.add(m_ann_opt.getStringValue());
        	}
        	command.add(m_annvcfdb.getStringValue());
        	command.add(vcf_infile);
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_annotated.vcf");
        	break;
		default:
			break;
        	
        }
        
        /**Intervals**/
        if(m_method.getStringValue().equals("Intervals")){
        	command.add("intervals");
        	command.add("-i "+vcf_infile);
        	if(m_interx.getBooleanValue()){
        		command.add("-x");
        	}
        	command.add(m_interbed.getStringValue());
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_intervals.vcf");
        }
        
        /**dbnsfp**/
        if(m_method.getStringValue().equals("Annotate with dbnsfp")){
        	command.add("dbnsfp");
        	if(m_dbnsfpfieldsall.getBooleanValue()){
        		command.add("-a");
        	}else{
        		if(m_dbnsfpfields.isActive() && m_dbnsfpfields.isEnabled()){
        			command.add("-f "+m_dbnsfpfields.getStringValue());
        		}	
        	}
        	command.add(m_dbnsfp.getStringValue());
        	command.add(vcf_infile);
        	stdOutFile = IO.replaceFileExtension(vcf_infile, ".snpSift_dbnfsp.vcf");
        }
        
    	/**Execute**/
        String [] cmd = new String [command.size()];
    	for(int i=0; i < cmd.length; i++) {
    		cmd[i] = command.get(i);
    	}
    	String lockFile = stdOutFile + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(cmd, exec, new File(lockFile),stdOutFile);
	
    	//Create Output Table
    	
    	if(stdOutFile.endsWith(".vcf")) {
    		outType = VCFCell.TYPE;
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL, outType).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(stdOutFile)};
    	
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

    	vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
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
    		if(CompatibilityChecker.inputFileNotOk(m_annvcfdb.getStringValue())) {
    			throw new InvalidSettingsException("Given annotation database invalid!");
    		}
    		outType = VCFCell.TYPE;
    		break;
    	default:
    		break;
    	}
    	//TODO check output according to method
    	
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL, outType).createSpec()})};
    }
}

