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



package de.helmholtz_muenchen.ibis.ngs.bwa;

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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;



/**
 * This is the model implementation of BWA.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiter
 */
public class BWANodeModel extends HTExecutorNodeModel {
    
	public static final String CFGKEY_REFSEQFILE 		= "refseqfile";
	public static final String CFGKEY_BWAFILE 	 		= "bwafile";
	public static final String CFGKEY_CHECKCOLORSPACED  = "checkColorSpaced";
	public static final String CFGKEY_BWTINDEX 			= "bwtIndex";
	public static final String CFGKEY_CHECKINDEX 		= "checkIndexRefSeq";
	public static final String CFGKEY_READGROUP 		= "readgroup";
	public static final String CFGKEY_READGROUPBOOLEAN  = "readgroupboolean";
	public static final String CFGKEY_ALNALGO 			= "alnalgo";
	public static final String CFGKEY_THREADS 			= "alnthreads";
	public static final String CFGKEY_OPTIONAL_Index 	= "optional_index";
	public static final String CFGKEY_OPTIONAL_Aln 		= "optional_aln";
	public static final String CFGKEY_OPTIONAL_Map 		= "optional_map";
	
	private final SettingsModelString m_refseqfile 				= new SettingsModelString(CFGKEY_REFSEQFILE,"");
	private final SettingsModelString m_bwafile 				= new SettingsModelString(CFGKEY_BWAFILE,"");
	private final SettingsModelBoolean m_checkIndexRefSeq 		= new SettingsModelBoolean(CFGKEY_CHECKINDEX,true);
//	private final SettingsModelBoolean m_checkColorSpaced		= new SettingsModelBoolean(CFGKEY_CHECKCOLORSPACED, false);
	private final SettingsModelString m_bwtIndex 				= new SettingsModelString(CFGKEY_BWTINDEX,"BWT-SW");
	private final SettingsModelString m_alnalgo 				= new SettingsModelString(CFGKEY_ALNALGO,"BWA-MEM");
	private final SettingsModelString m_readGroup 				= new SettingsModelString(CFGKEY_READGROUP,"@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA");
	private final SettingsModelBoolean m_readGroupBoolean 		= new SettingsModelBoolean(CFGKEY_READGROUPBOOLEAN,false);
	private final SettingsModelIntegerBounded m_ALN_THREADS 	= new SettingsModelIntegerBounded(CFGKEY_THREADS,2, 1, Integer.MAX_VALUE);
//	private final SettingsModelOptionalString m_Optional_Index 	= new SettingsModelOptionalString(CFGKEY_OPTIONAL_Index,"",false);
	private final SettingsModelOptionalString m_Optional_Aln 	= new SettingsModelOptionalString(CFGKEY_OPTIONAL_Aln,"",false);
	private final SettingsModelOptionalString m_Optional_Map 	= new SettingsModelOptionalString(CFGKEY_OPTIONAL_Map,"",false);
	

	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BWANodeModel.class);
	private static String readType = "";
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SAMFile";
	
	
    /**
     * Constructor for the node model.
     */
    protected BWANodeModel() {
    	
    	super(1, 1);
    	addSetting(m_bwafile);
    	addSetting(m_refseqfile);
    	addSetting(m_bwtIndex);
//    	addSetting(m_checkColorSpaced);
    	addSetting(m_checkIndexRefSeq);
    	addSetting(m_readGroup);
    	addSetting(m_readGroupBoolean);
    	addSetting(m_alnalgo);
    	addSetting(m_ALN_THREADS);
//    	addSetting(m_Optional_Index);
    	addSetting(m_Optional_Aln);
    	addSetting(m_Optional_Map);
    	
    }

   
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	
    	/**
		 * Get the Parameters
		 */
    	String path2refFile 	= m_refseqfile.getStringValue();
    	String path2readFile 	= inData[0].iterator().next().getCell(0).toString();
    	String path2readFile2 	= "";
 
    	if(readType.equals("paired-end")){
    		path2readFile2 = inData[0].iterator().next().getCell(1).toString();	
    	}
 
    	String basePath 		= path2readFile.substring(0,path2readFile.lastIndexOf('/')+1);
    	String outBaseName1 	= path2readFile.substring(path2readFile.lastIndexOf("/")+1,path2readFile.lastIndexOf("."));
    	String outBaseName 		= outBaseName1;
    	String outBaseName2 	= outBaseName1;
    	String memOut 			= basePath+outBaseName1+"_mem.sam";
    	String path2bwa 		= m_bwafile.getStringValue();
    	int threads 			= m_ALN_THREADS.getIntValue();
    	
    	if(readType.equals("paired-end")){
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2readFile.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}
    	String out2Name  = basePath+outBaseName+"_aln.sam";
    	String out1Name  = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";

    	Boolean isBam = false;	
    	if(path2readFile.substring(path2readFile.length()-3, path2readFile.length()) == "bam") {
    		path2readFile2 	= path2readFile;
    		isBam 			= true;
    	}	

    	//Prepare Index
    	bwa_index(exec,path2bwa, path2refFile);

    	//BWA aln
    	if(m_alnalgo.getStringValue().equals("BWA-backtrack")){
        	LOGGER.info("Find the SA coordinates of the input reads.\n");
        	bwa_aln(exec,readType, basePath, outBaseName, outBaseName1, outBaseName2, path2refFile, path2bwa, path2readFile, path2readFile2, isBam,threads);
        	LOGGER.info("Finished BWA aln...");
    	}
    	//BWA Mapping
    	bwa_map(exec, readType,path2bwa,path2refFile,path2readFile,out1Name,out2Name,out11Name,out12Name,path2readFile2,memOut,threads);
    	
     	
    	/**
    	 * OUTPUT
    	 */
    	if(m_alnalgo.getStringValue().equals("BWA-MEM")){
    		out2Name = memOut;
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    			new DataColumnSpecCreator(OUT_COL1, SAMCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{FileCellFactory.create(out2Name)};
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	   	
    	
		return new BufferedDataTable[]{outTable};
    }
    

    /**
     * Runs bwa index
     * @param lockFile
     * @param colorSpaced
     * @param path2bwa
     * @param path2refFile
     * @param path2readFile
     * @throws Exception 
     */
    private void bwa_index(ExecutionContext exec,String path2bwa, String path2refFile) throws Exception{
    	/**Only execute if Index needs to be created**/
    	if(m_checkIndexRefSeq.getBooleanValue()) {
    		
    		LOGGER.info("Indexing reference sequence.\n");
    		
        	ArrayList<String> command = new ArrayList<String>();
        	// Constant values
        	command.add(path2bwa+" index");
        	      	
        	//Indexing Type
	    	if(m_bwtIndex.getStringValue().equals("BWT-SW")) {
	    		command.add("-a bwtsw");
	    	} else if(m_bwtIndex.getStringValue().equals("IS")){
	    		command.add("-a is");
	    	}else{
	    		throw new InvalidSettingsException("Oh oh oh! No valid indexing algorithm!");
	    	}
//	    	// Colorspace
//	    	if(m_checkColorSpaced.getBooleanValue()) {
//	    		command.add("-c");
//	    	}
	    	
	    	//Add Reference genome
	    	command.add(path2refFile);

	    	/**Execute**/
	    	String lockFile = path2refFile + SuccessfulRunChecker.LOCK_ENDING;
	    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));
			
    	} else {
    		LOGGER.info("Indexing reference sequence SKIPPED.\n");
    	}
    }
    
    private void bwa_aln(ExecutionContext exec, String readType, String basePath, String outBaseName, String outBaseName1, String outBaseName2, String path2refFile, String path2bwa, String path2readFile, String path2readFile2, boolean isBam, int Threads) throws Exception{
    	
    	ArrayList<String> command = new ArrayList<String>();
    	// Constant values
    	command.add(path2bwa+" aln");
    	
    	String outName = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";
    	String outfile = outName;
    	
    	//Multi-Threading 
    	command.add("-t " +Threads);
    	
    	command.add(m_Optional_Aln.getStringValue());
    	
    	//If Inputfile is in bam format
    	if(readType.equals("paired-end")){
    		outfile = out11Name;				//Set Outfile for forward reads
    		if(isBam){
    			command.add("-b1");
    		}
    	}else{	// Single-end
    		if(isBam){
    			command.add("-b0");
    		}
		}
    	
    	//Perform aln for forward reads OR single end reads
    	command.add(path2refFile);
    	command.add(path2readFile);
    	command.add("-f "+outfile);
    	
    	String lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    	/**Execute**/
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));
   
		//If paired end, repeat previous step
    	if(readType.equals("paired-end")) {
    		
        	if(isBam){
        		command.set(2, "-b2");        	
        		command.set(4, path2readFile2);
            	command.set(5, " -f "+ out12Name);
        	}else{
        		command.set(3, path2readFile2);
            	command.set(4, " -f "+ out12Name);
        	}
        	/**Execute**/
        	lockFile = out12Name + SuccessfulRunChecker.LOCK_ENDING;
        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));
		}
    }
    
    
  private void bwa_map(ExecutionContext exec, String readType, String path2bwa, String path2refFile, String path2readFile, String out1Name, String out2Name, String out11Name, String out12Name, String path2readFile2, String memOut, int threads) throws Exception{ 	
    	
  		ArrayList<String> command = new ArrayList<String>();
  		String alnalgo = m_alnalgo.getStringValue(); 	
    	
  		// ### BWA-Backtrack ###
    	if(alnalgo.equals("BWA-backtrack")) {
    		
        	if(readType.equals("single-end")) {
        		// bwa samse sequence.fasta aln_sa.sai s_1_1_sequence.txt > aln.sam
        		LOGGER.info("Generate alignments in the SAM format given single-end reads.\n");
        		command.add(path2bwa+" samse");
        		
        		command.add(m_Optional_Map.getStringValue());
        		
        		/**Other Options**/
        		if(m_readGroupBoolean.getBooleanValue()){
        			command.add("-r "+m_readGroup.getStringValue());
        		}
        		/**In and Outfiles**/
        		command.add("-f "+out2Name);
        		command.add(path2refFile);
        		command.add(out1Name);
        		command.add(path2readFile);

        	} else {
        		// bwa sampe sequence.fasta aln_sa_1.sai aln_sa_2.sai s_1_1_sequence.fq s_1_2_sequence.fq > aln.sam
        		LOGGER.info("Generate alignments in the SAM format given paired-end reads.\n");
        		command.add(path2bwa+" sampe");
        		
        		command.add(m_Optional_Map.getStringValue());
        		
        		/**Other Options**/
        		if(m_readGroupBoolean.getBooleanValue()){
        			command.add("-r "+m_readGroup.getStringValue());
        		}
        		/**In and Outfiles**/
        		command.add("-f "+out2Name);
        		command.add(path2refFile);
        		command.add(out11Name);
        		command.add(out12Name);
        		command.add(path2readFile);
        		command.add(path2readFile2);
        		
        	}
// ### BWA-SW ###
    	} else if(alnalgo.equals("BWA-SW")) {
        	LOGGER.info("Generate alignments in the SAM format.\n");

        	command.add(path2bwa+" bwasw");
        	command.add("-t "+threads);
        	
        	command.add(m_Optional_Map.getStringValue());
        	
    		/**In and Outfiles**/
        	command.add("-f "+out2Name);
        	command.add(path2refFile);
        	command.add(path2readFile);
    		if(readType.equals("paired-end")) {
    			command.add(path2readFile2);
    		} 		
// ### BWA-MEM ###
    	} else if(alnalgo.equals("BWA-MEM")) {
        	LOGGER.info("Generate alignments in the SAM format.\n");

        	command.add(path2bwa+" mem");
        	command.add("-t "+threads);
        	command.add(m_Optional_Map.getStringValue());
        	
    		if(m_readGroupBoolean.getBooleanValue()){
    			command.add("-R "+m_readGroup.getStringValue());
    		}
        	
    		/**In and Outfiles**/
        	command.add(path2refFile);
        	command.add(path2readFile);
    		if(readType.equals("paired-end")) {
    			command.add(path2readFile2);
    		} 
    	}
		
    	String lockFile = out2Name + SuccessfulRunChecker.LOCK_ENDING;
		/** Execute **/
		if (alnalgo.equals("BWA-MEM")) {
			String stdErr = IO.replaceFileExtension(memOut,"stdErr");
			
			super.executeCommand(new String[] { StringUtils.join(command, " ") }, exec, null, new File(lockFile), memOut,stdErr, null, null, null);
		} else {
			super.executeCommand(new String[] { StringUtils.join(command, " ") }, exec,new File(lockFile));
		}
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
    	
    	
		if(CompatibilityChecker.inputFileNotOk(m_bwafile.getStringValue(), false)) {
			throw new InvalidSettingsException("Set path to BWA binary!");
		}
    	
    	    	
        //Version control
    	try {
	        if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==1){
	        	setWarningMessage("WARNING: You are using a newer BWA version than "+FileValidator.BWA_VERSION +"! This may cause problems!");
	        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==2){
	        	setWarningMessage("WARNING: You are using an older BWA version than "+FileValidator.BWA_VERSION +"! This may cause problems!");
	        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==-1){
	        	setWarningMessage("Your BWA version could not be determined! Correct behaviour can only be ensured for BWA version "+FileValidator.BWA_VERSION+".");
	        }
    	} catch (Exception e) {
    		throw new InvalidSettingsException("Specify a valid BWA version!");
    	}
   
    	
    	if(m_refseqfile.getStringValue().length() > 1) {
	    	if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	            throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
	    	}
    	}
    	

        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, SAMCell.TYPE).createSpec()})};
    }
}

