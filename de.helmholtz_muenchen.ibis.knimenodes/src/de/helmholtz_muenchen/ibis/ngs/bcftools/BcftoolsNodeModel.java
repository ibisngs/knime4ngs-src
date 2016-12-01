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
package de.helmholtz_muenchen.ibis.ngs.bcftools;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
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
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

/**
 * This is the model implementation of Bcftools.
 * 
 *
 * @author Maximilian Hastreiter

 */
public class BcftoolsNodeModel extends HTExecutorNodeModel {

	
	/**
	 * Bcftools method chooser and infile
	 */
	public static final String CFGKEY_BCFMETHOD="bcfmethod";
	public static final String CFGKEY_PATH2BCFTOOLS="path2bcftools";
	public static final String CFGKEY_VCFSAMPLEHEADER="VCFSAMPLEHEADER";
	
	private final SettingsModelString m_path2bcftools = new SettingsModelString(
			CFGKEY_PATH2BCFTOOLS, "");
	private final SettingsModelString m_bcfmethod = new SettingsModelString(
			CFGKEY_BCFMETHOD, "");
	
	
	/**
	 * Concat
	 */
	public static final String CFGKEY_CONCAT_OVERLAPS 		= "concat_overlaps";
	public static final String CFGKEY_CONCAT_OUTFILE_TYPE 	= "concat_outfile_type";
	private final SettingsModelBoolean m_concat_overlap 	= new SettingsModelBoolean(
			BcftoolsNodeModel.CFGKEY_CONCAT_OVERLAPS, false);
	private final SettingsModelString m_concat_outfile_type = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_CONCAT_OUTFILE_TYPE, "");
	
	private String bcf_bin;
	
	/**
	 * Call
	 */
//	public static final String CFGKEY_CALL_OUTFILE_TYPE 	= "call_outfile_type";
//	public static final String CFGKEY_KEEPALLELES           = "keepalleles";
//	public static final String CFGKEY_BEDFILE				= "bedfile";
//	public static final String CFGKEY_SAMPLELIST			= "samplelist";
//	public static final String CFGKEY_IFBEDFILE				= "ifbedfile";
//	public static final String CFGKEY_IFSAMPLELIST			= "ifsamplelist";
//	public static final String CFGKEY_SNPCALLING			= "snpcalling";
//	public static final String CFGKEY_CONSTRAINEDCALLING	= "constrainedcalling";
//	public static final String CFGKEY_OUTVARIANTSONLY		= "outvariantsonly";
//	
//	
//	private final SettingsModelString m_call_outfile_type = new SettingsModelString(
//			BcftoolsNodeModel.CFGKEY_CALL_OUTFILE_TYPE, "");
//	private final SettingsModelBoolean m_keepallel = new SettingsModelBoolean(
//			BcftoolsNodeModel.CFGKEY_KEEPALLELES, false);
//	private final SettingsModelString m_bedfile = new SettingsModelString(
//			BcftoolsNodeModel.CFGKEY_BEDFILE,"");
//	private final SettingsModelBoolean m_ifsamplelist = new SettingsModelBoolean(
//			BcftoolsNodeModel.CFGKEY_IFSAMPLELIST, false);
//	private final SettingsModelString m_samplelist = new SettingsModelString(
//			BcftoolsNodeModel.CFGKEY_SAMPLELIST,"");
//	private final SettingsModelBoolean m_ifbedfile = new SettingsModelBoolean(
//			CFGKEY_IFBEDFILE, false);
//	private final SettingsModelBoolean m_snpcalling = new SettingsModelBoolean(
//			CFGKEY_SNPCALLING, false);
//	private final SettingsModelString m_constrainedcalling = new SettingsModelString(
//			BcftoolsNodeModel.CFGKEY_CONSTRAINEDCALLING,"");
//	private final SettingsModelBoolean m_outvariantsonly = new SettingsModelBoolean(
//			CFGKEY_OUTVARIANTSONLY, false);
	
	
	
	//Reheader
	private final SettingsModelString m_vcfsampleheader = new SettingsModelString(
			CFGKEY_VCFSAMPLEHEADER, "");
	
	//Extra
	public static final String CFGKEY_FURTHER_OPTIONS 	= "further_options";
	private final SettingsModelOptionalString m_furtherOptions = new SettingsModelOptionalString(
	BcftoolsNodeModel.CFGKEY_FURTHER_OPTIONS,"",false);
	
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(BcftoolsNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2Bcftools";
	public static final String OUT_COL2 = "Path2BcftoolsOutput";
	
	protected String OUTFILE;
	
    /**
     * Constructor for the node model.
     */
    protected BcftoolsNodeModel() {
    
        super(1,1);
        
        addSetting(m_bcfmethod);
   	 	addSetting(m_path2bcftools);
   	 	addSetting(m_vcfsampleheader);
   	 	addSetting(m_furtherOptions);
	   	//Concat
	   	addSetting(m_concat_outfile_type);
	    addSetting(m_concat_overlap);
	    
	    addPrefPageSetting(m_path2bcftools, IBISKNIMENodesPlugin.BCFTOOLS);
        
//        m_bedfile.setEnabled(false);
//        m_samplelist.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	ArrayList<String> command = new ArrayList<String>();
    	String method = m_bcfmethod.getStringValue();
    	command.add(bcf_bin+" "+method);

//    	if(m_bcfmethod.getStringValue().equals("call")){
//    		command.addAll(call(inData));
//    	}
    	if(m_bcfmethod.getStringValue().equals("concat")){
    		command.addAll(concat(inData));
    	}
    	
    	else if(m_bcfmethod.getStringValue().equals("index")){
            command.addAll(index(inData));
    	}
    	
    	else if(m_bcfmethod.getStringValue().equals("reheader")){
    		 command.addAll(reheader(inData));
    	}
    	
    	else if(m_bcfmethod.getStringValue().equals("stats")){
    		command.addAll(stats(inData));
    	}
   	
    	
    	/**
    	 * Execute
    	 */
    	
    	String lockFile = OUTFILE +"_"+m_bcfmethod.getStringValue()+SuccessfulRunChecker.LOCK_ENDING;
		super.executeCommand(new String[]{StringUtils.join(command, " ")}, OUTFILE ,exec,new File(lockFile),OUTFILE,OUTFILE+".stdErr");
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{FileCellFactory.create(OUTFILE)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
   
		return new BufferedDataTable[]{outTable};
    }
    
    /**
     * bcftools call
     * @param inData
     * @return Command options for bcftools call
     * @throws InvalidSettingsException
     */
//    protected ArrayList<String> call(BufferedDataTable[] inData) throws InvalidSettingsException{
//    	//Row Iterator
//    	Iterator<DataRow> it = inData[0].iterator();
//    	ArrayList<String> inputfiles = new ArrayList<String>();
//        // Get input files
//        while(it.hasNext()){
//        	DataRow curr_Row = it.next();
//        	
//        	//File from current row
//        	String curr_File = curr_Row.getCell(0).toString();
//        	//Check format and append if correct
//        	if(curr_File.endsWith("vcf.gz") || curr_File.endsWith("bcf")){
//        		
//        		if(!(new File(curr_File+".tbi").exists()) && !curr_File.endsWith("bcf")){
//        			throw new InvalidSettingsException("Missing tabix index for "+curr_File);
//        		}
//        		inputfiles.add(curr_File);
//        	}else{
//        		throw new InvalidSettingsException("Cannot work with current input file. Required format: vcf.gz | bcf.gz");
//        	}
//        	
//        }
//    	
//       
//        ArrayList<String> command = new ArrayList<String>();
//   	 
//    	if(m_keepallel.getBooleanValue()){
//    		command.add("-A");
//    	}
//   	
//    	if(m_bedfile.isEnabled()){
//    		command.add("-R "+m_bedfile.getStringValue());
//    	}
//
//    	if(m_samplelist.isEnabled()){
//    		command.add("-s "+m_samplelist.getStringValue());
//    	}
//	
//    	String outfile = inputfiles.get(0).split(".vcf")[0]+".call.";
//        String outfile_type = m_call_outfile_type.getStringValue();
//        
//        if(outfile_type.equals("compressed BCF")){
//        	command.add("--output-type b");
//        	outfile+="bcf.gz";
//        }else if(outfile_type.equals("uncompressed BCF")){
//        	command.add("--output-type u");
//        	outfile+="bcf";
//        }else if(outfile_type.equals("compressed VCF")){
//        	command.add("--output-type z");
//        	outfile+="vcf.gz";
//        }else if(outfile_type.equals("uncompressed VCF")){
//        	command.add("--output-type v");
//        	outfile+="vcf";
//        }
//        command.add("--output "+outfile);
//        OUTFILE=outfile;
//
//        if(m_snpcalling.getBooleanValue()){
//        	command.add("-c");
//        }
//        
//        if(!m_constrainedcalling.getStringValue().equals("No constrains")){
//        	command.add("-C "+m_constrainedcalling.getStringValue());
//        }
//        if(m_outvariantsonly.getBooleanValue()){
//        	command.add("-v");
//        }
//
		//Add extra Options
//		if(m_furtherOptions.isEnabled()){
//			command.add(m_furtherOptions.getStringValue());
//		}
//    
//        command.add(inputfiles.get(0));
//        return command;
//    }

    /**
     * bcftools reheader
     * @param inData
     * @return Command options for bcftools reheader
     * @throws InvalidSettingsException
     */
    protected ArrayList<String> reheader(BufferedDataTable[] inData) throws InvalidSettingsException{
    	//Row Iterator
    	Iterator<DataRow> it = inData[0].iterator();

    	ArrayList<String> inputfiles = new ArrayList<String>();
        // Get input files
        while(it.hasNext()){
        	DataRow curr_Row = it.next();
        	
        	//File from current row
        	String curr_File = curr_Row.getCell(0).toString();
        	//Check format and append if correct
        	if(curr_File.endsWith("vcf.gz") || curr_File.endsWith("bcf.gz") || curr_File.endsWith("vcf") || curr_File.endsWith("bcf")){
        		
//        		if(!new File(curr_File+".tbi").exists()){
//        			throw new InvalidSettingsException("Missing tabix index for "+curr_File);
//        		}
        		inputfiles.add(curr_File);
        	}else{
        		throw new InvalidSettingsException("Cannot work with current input files. Required format: vcf.gz");
        	}
        	
        }
    	
        ArrayList<String> command = new ArrayList<String>();
        
        command.add("--samples "+m_vcfsampleheader.getStringValue());
        
        String outfile = inputfiles.get(0);  
        
        if(inputfiles.get(0).endsWith(".vcf")){
        	outfile = outfile.replace(".vcf", "_reheader.vcf");
        }else if(inputfiles.get(0).endsWith(".vcf.gz")){
        	outfile = outfile.replace(".vcf.gz", "_reheader.vcf.gz");
        }else if(inputfiles.get(0).endsWith(".bcf.gz")){
        	outfile = outfile.replace(".bcf.gz", "_reheader.bcf.gz");
        }else if(inputfiles.get(0).endsWith(".bcf")){
        	outfile = outfile.replace(".bcf", "_reheader.bcf");
        }
      
        command.add("--output "+outfile);
        OUTFILE = outfile;
        command.add(inputfiles.get(0));
        
        return command;
    }
    
    
    
    /**
     * bcftools index
     * @param inData
     * @return Command options for bcftools index
     * @throws InvalidSettingsException
     */
    protected ArrayList<String> stats(BufferedDataTable[] inData) throws InvalidSettingsException{
    	//Row Iterator
    	Iterator<DataRow> it = inData[0].iterator();

    	ArrayList<String> inputfiles = new ArrayList<String>();
        // Get input files
        while(it.hasNext()){
        	DataRow curr_Row = it.next();
        	
        	//File from current row
        	String curr_File = curr_Row.getCell(0).toString();
        	//Check format and append if correct
        	if(curr_File.endsWith("vcf.gz")){
        		
        		if(!new File(curr_File+".tbi").exists()){
        			throw new InvalidSettingsException("Missing tabix index for "+curr_File);
        		}
        		inputfiles.add(curr_File);
        	}else{
        		throw new InvalidSettingsException("Cannot work with current input files. Required format: vcf.gz");
        	}
        	
        }
    	
        ArrayList<String> command = new ArrayList<String>();
        
    	//Add extra Options
    	if(m_furtherOptions.isActive()){
    		command.add(m_furtherOptions.getStringValue());
    	}
        
        command.add(inputfiles.get(0));
        if(inputfiles.size()>1){
        	command.add(inputfiles.get(1));
        }
        
        OUTFILE=inputfiles.get(0)+".vchk";
        
        return command;
    }
    
    
    /**
     * bcftools index
     * @param inData
     * @return Command options for bcftools index
     * @throws InvalidSettingsException
     */
    protected ArrayList<String> index(BufferedDataTable[] inData) throws InvalidSettingsException{
    	//Row Iterator
    	Iterator<DataRow> it = inData[0].iterator();

    	ArrayList<String> inputfiles = new ArrayList<String>();
        // Get input files
        while(it.hasNext()){
        	DataRow curr_Row = it.next();
        	
        	//File from current row
        	String curr_File = curr_Row.getCell(0).toString();
        	//Check format and append if correct
        	if(curr_File.endsWith("vcf.gz") || curr_File.endsWith("bcf.gz")){
        		
        		if(!new File(curr_File+".tbi").exists()){
        			throw new InvalidSettingsException("Missing tabix index for "+curr_File);
        		}
        		inputfiles.add(curr_File);
        	}else{
        		throw new InvalidSettingsException("Cannot work with current input files. Required format: vcf.gz | bcf.gz");
        	}
        	
        }
    	
        ArrayList<String> command = new ArrayList<String>();
    	//Add extra Options
    	if(m_furtherOptions.isActive()){
    		command.add(m_furtherOptions.getStringValue());
    	}
        command.add(inputfiles.get(0));     
        OUTFILE=inputfiles.get(0)+".csi";
        
        return command;
    }

    /**
     * Creates command for bcftools concat
     * @param inData
     * @return Command options for bcftools concat
     * @throws InvalidSettingsException
     */
    protected ArrayList<String> concat(BufferedDataTable[] inData) throws InvalidSettingsException{
    
    	//Row Iterator
    	Iterator<DataRow> it = inData[0].iterator();

    	ArrayList<String> inputfiles = new ArrayList<String>();
        // Get input files
        while(it.hasNext()){
        	DataRow curr_Row = it.next();
        	
        	//File from current row
        	String curr_File = curr_Row.getCell(0).toString();
        	//Check format and append if correct
        	if(curr_File.endsWith("vcf.gz")){
        		
        		if(!new File(curr_File+".tbi").exists()){
        			throw new InvalidSettingsException("Missing tabix index for "+curr_File);
        		}
        		inputfiles.add(curr_File);
        	}else{
        		throw new InvalidSettingsException("Cannot work with current input files. Required format: vcf.gz");
        	}
        	
        }
    	
        ArrayList<String> command = new ArrayList<String>();
        
        if(m_concat_overlap.getBooleanValue()){
        	command.add("--allow-overlaps");
        }
        
        String outfile = inputfiles.get(0).split(".vcf")[0]+".concat.";
        String outfile_type = m_concat_outfile_type.getStringValue();
        
        if(outfile_type.equals("compressed BCF")){
        	command.add("--output-type b");
        	outfile+="bcf.gz";
        }else if(outfile_type.equals("uncompressed BCF")){
        	command.add("--output-type u");
        	outfile+="bcf";
        }else if(outfile_type.equals("compressed VCF")){
        	command.add("--output-type z");
        	outfile+="vcf.gz";
        }else if(outfile_type.equals("uncompressed VCF")){
        	command.add("--output-type v");
        	outfile+="vcf";
        }
        
    	//Add extra Options
    	if(m_furtherOptions.isActive()){
    		command.add(m_furtherOptions.getStringValue());
    	}
        
        command.add("--output "+outfile);
        command.addAll(inputfiles);
        OUTFILE=outfile;
        
    	return command;
    }
    
    
    
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

		super.updatePrefs();
		bcf_bin = m_path2bcftools.getStringValue();
		
		if (CompatibilityChecker.inputFileNotOk(bcf_bin, false)) {
			throw new InvalidSettingsException("Invalid path to bcftools!");
		}

		return new DataTableSpec[] { new DataTableSpec(
				new DataColumnSpec[] { new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec() }) };

	}
}

