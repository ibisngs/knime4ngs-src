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
package de.helmholtz_muenchen.ibis.ngs.samtools;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;


/**
 * This is the model implementation of SamTools.
 * 
 * @author Sebastian Kopetzky
 * @author Maximilian Hastreiter
 */

public class SamToolsNodeModel extends HTExecutorNodeModel {
    
	
	protected static boolean useSamtools = true;
	protected static boolean useBamfile = true;
	protected static boolean useFastafile = true;
	
	public static final String CFGKEY_UTILITY = "utility";
	public static final String CFGKEY_SAMTOOLS = "samtools";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_OPTIONAL = "optional";
	
	
	private final SettingsModelString m_utility 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_UTILITY, "");
	private final SettingsModelString m_samtools 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_SAMTOOLS, "");
	private final SettingsModelString m_refseqfile 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_REFSEQFILE, "");
	private final SettingsModelOptionalString m_Optional 	= new SettingsModelOptionalString(CFGKEY_OPTIONAL,"",false);


	/**Parameters for "calmd"**/
//	public static final String CFGKEY_CHANGEIDENTBASES = "changeIdentBases";
//	public static final String CFGKEY_USECOMPRESSION = "useCompression";
//	public static final String CFGKEY_COMPRESSION = "compression";
//	public static final String CFGKEY_INPUTISSAM = "inputIsSam";
//	public static final String CFGKEY_MODIFYQUAL = "modifyQual";
//	public static final String CFGKEY_BQTAG = "bqTag";
//	public static final String CFGKEY_EXTENDEDBAQ = "extendedBAQ";
//	public static final String CFGKEY_DOCAPMAPQUAL = "doCapMapQual";
//	public static final String CFGKEY_CAPMAPQUAL = "capMapQual";
	
//	private final SettingsModelBoolean m_changeIdentBases = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_CHANGEIDENTBASES, false);
//	private final SettingsModelBoolean m_useCompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USECOMPRESSION, false);	
//	private final SettingsModelString m_compression = new SettingsModelString(SamToolsNodeModel.CFGKEY_COMPRESSION, "");
//	private final SettingsModelBoolean m_inputIsSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_INPUTISSAM, false);
//	private final SettingsModelBoolean m_modifyQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MODIFYQUAL, false);
//	private final SettingsModelBoolean m_bqTag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_BQTAG, false);
//	private final SettingsModelBoolean m_extendedBAQ = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_EXTENDEDBAQ, false);

	/**Parameters for rmdup**/
	public static final String CFGKEY_REMOVEDUP = "removeDup";
	public static final String CFGKEY_TREATPE = "treatPE";
	private final SettingsModelBoolean m_removeDup = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_REMOVEDUP, false);
	private final SettingsModelBoolean m_treatPE = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_TREATPE, false);
	/**Parameters for cat**/
	public static final String CFGKEY_USEHEADERSAM = "useHeaderSAM";
	public static final String CFGKEY_HEADERSAM = "headerSAM";

	private final SettingsModelBoolean m_useHeaderSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEHEADERSAM, false);
	private final SettingsModelString m_headerSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_HEADERSAM, "");

	/**Parameters for reheader**/
	public static final String CFGKEY_REHINSAM = "rehInSAM";

	private final SettingsModelString m_rehInSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINSAM, "");
	
	/**Parameters for merge**/
	public static final String CFGKEY_MCOMPRESSION = "mcompression"; //zlib compression
	public static final String CFGKEY_MFORCE = "mforce";
	public static final String CFGKEY_USEMHFILE = "usemhfile";
	public static final String CFGKEY_MHFILE = "mhfile";
	public static final String CFGKEY_MSORTED = "msorted";
	public static final String CFGKEY_MREGION = "mregion";
	public static final String CFGKEY_USEMREGION = "usemregion";
	public static final String CFGKEY_MRGTAG = "mrgtag";
	public static final String CFGKEY_MUNCOMPRESSED = "muncompressed";
	
	private final SettingsModelBoolean m_mcompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MCOMPRESSION, false);
	private final SettingsModelBoolean m_mforce = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MFORCE, false);
	private final SettingsModelBoolean m_usemhfile = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMHFILE, false);
	private final SettingsModelString m_mhfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_MHFILE, "");
	private final SettingsModelBoolean m_msorted = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MSORTED, false);
	private final SettingsModelString m_mregion = new SettingsModelString(SamToolsNodeModel.CFGKEY_MREGION, "");
	private final SettingsModelBoolean m_usemregion = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMREGION, false);
	private final SettingsModelBoolean m_mrgtag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MRGTAG, false);
	private final SettingsModelBoolean m_muncompressed = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MUNCOMPRESSED, false);
	/**faidx**/

	/**phase**/
	public static final String CFGKEY_BLOCKLENGTH = "blocklength";
	public static final String CFGKEY_PREFIX = "prefix";
	public static final String CFGKEY_HETPHRED = "hetphred";
	public static final String CFGKEY_MINQUAL = "minqual";
	public static final String CFGKEY_MAXDEPTH = "maxdepth";
	public static final String CFGKEY_FIXCHIMERAS = "fixchimeras";

	private final SettingsModelIntegerBounded m_blocklength = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_BLOCKLENGTH, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelString m_prefix = new SettingsModelString(SamToolsNodeModel.CFGKEY_PREFIX,"phase");
	private final SettingsModelIntegerBounded m_hetphred = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_HETPHRED, 37, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minqual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MINQUAL, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxdepth = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MAXDEPTH, 256, 1, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_fixchimeras = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_FIXCHIMERAS, false);

	public static String OUT_COL1 			= "Outfile";
	private String OutCellType 				= "FileCell";	
	private String samtools_bin, ref_genome;
	
	/**
     * Constructor for the node model.
     */
    protected SamToolsNodeModel() {
    
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1), 1);
       
       	addSetting(m_utility);
    	addSetting(m_samtools);
//    	addSetting(m_bamfile);
    	addSetting(m_refseqfile);
    	addSetting(m_Optional);
    	
    	//calmd
//    	addSetting(m_changeIdentBases);
//    	addSetting(m_compression);
//    	addSetting(m_useCompression);
//    	addSetting(m_inputIsSAM);
//    	addSetting(m_modifyQual); 
//    	addSetting(m_bqTag); 
//    	addSetting(m_extendedBAQ); 

    	//remdup
    	addSetting(m_removeDup);
    	addSetting(m_treatPE);
    	//cat
    	addSetting(m_useHeaderSAM);
    	addSetting(m_headerSAM);
//    	addSetting(m_inBAM1);

    	//reheader
    	addSetting(m_rehInSAM);

    	//merge
    	addSetting(m_mcompression);
    	addSetting(m_mforce);
    	addSetting(m_usemhfile);
    	addSetting(m_mhfile);
    	addSetting(m_msorted);
    	addSetting(m_mregion);
    	addSetting(m_usemregion);
    	addSetting(m_mrgtag);
    	addSetting(m_muncompressed);

    	//phase
    	addSetting(m_blocklength);
    	addSetting(m_prefix);
    	addSetting(m_hetphred);
    	addSetting(m_minqual);
    	addSetting(m_maxdepth);
    	addSetting(m_fixchimeras);
    	
		addPrefPageSetting(m_samtools, IBISKNIMENodesPlugin.SAMTOOLS);
		addPrefPageSetting(m_refseqfile, IBISKNIMENodesPlugin.REF_GENOME);

    }

    /**
     * {@inheritDoc}
     * @return 
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	String lockFile 		= "";
    	String OutCellType 		= "FileCell";	
    	
    	int colIndex = 0;
    	if(SamToolsNodeDialog.getUseMainInputColBool()){
    		colIndex = inData[0].getDataTableSpec().findColumnIndex(SamToolsNodeDialog.getMainInputCol1());
    	}
    
    	String path2bamfile = inData[0].iterator().next().getCell(colIndex).toString();
    	
    	String basePathWithFileName = "";
    	if(!path2bamfile.isEmpty()){
    	   	basePathWithFileName = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));	
    	}
    	else{
    		basePathWithFileName = ref_genome.substring(0,ref_genome.lastIndexOf("."));	
    	}
 
    	String utility = m_utility.getStringValue();
    	ArrayList<String> command = new ArrayList<String>();
    	command.add(samtools_bin+" "+utility);
    	
    	if(m_Optional.isActive()){
    		command.add(m_Optional.getStringValue());	
    	}
    	
    	String outfile = "-1";
    	boolean stdOut = false;
    	
    	
    	if(utility.equals("flagstat") || utility.equals("idxstats") || utility.equals("depth")) {
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".txt";
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		stdOut = true;
    		
//    	} else if(utility.equals("calmd")) {
//
//    		String outFileExtension = ".sam";
//    		OutCellType = "SAMCell";
//    		
//    		if(m_changeIdentBases.getBooleanValue()){
//    			command.add("-e");
//    		}
//    		if(m_useCompression.getBooleanValue()){
//    			outFileExtension = ".bam";
//    			OutCellType = "BAMCell";
//    			if(m_compression.getStringValue().equals("uncompressed")){
//    				command.add("-u");
//    			}
//    			else{
//    				command.add("-b");
//    			}
//    		}
//    		if(m_inputIsSAM.getBooleanValue()){
//    			command.add("-S");
//    		}
//    		if(m_modifyQual.getBooleanValue()){
//    			command.add("-A");
//    		}
//    		if(m_bqTag.getBooleanValue()){
//    			command.add("-r");
//    		}
//    		if(m_extendedBAQ.getBooleanValue()){
//    			command.add("-E");
//    		}
//    		
//    		command.add(path2bamfile);
//    		command.add(path2seqfile);
//    		outfile = basePathWithFileName + "_" + utility + outFileExtension;
//    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
//    		stdOut = true;
    		
    		
    			
    	} else if (utility.equals("fixmate")) {
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".bam";
    		command.add(outfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("rmdup")) {
    		if(m_removeDup.getBooleanValue()){
    			command.add("-s");
    		}
    		if(m_treatPE.getBooleanValue()){
    			command.add("-S");
    		}
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".bam";
    		command.add(outfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("cat")) {
    		
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += ".catAll.bam";
    		if(m_useHeaderSAM.getBooleanValue()){
    			command.add("-h " + IO.processFilePath(m_headerSAM.getStringValue()));
    		}
    		command.add("-o " + outfile);
    		
            Iterator<DataRow> it = inData[0].iterator();
            int size = 0;
            while(it.hasNext()){
            	command.add(it.next().getCell(colIndex).toString());
            	size++;
            } 
     		if(size<2){
     			throw new InvalidSettingsException("Insufficient number of input files. At least 2 bam files required for samtools cat.");
     		}
     		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
     		OutCellType = "BAMCell";
     		
     		
    	} else if (utility.equals("faidx")){
    		command.add(ref_genome);
    		outfile = ref_genome;
    		lockFile = outfile+".faidx" + SuccessfulRunChecker.LOCK_ENDING;
    		
    	} else if (utility.equals("reheader")) {
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")) + "_" + utility + ".bam";
    		command.add(IO.processFilePath(m_rehInSAM.getStringValue()));
    		command.add(path2bamfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		stdOut = true;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("merge")) {
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += ".mergedAll.bam";
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		
    		if(m_mcompression.getBooleanValue()){
    			command.add("-1");    			
    		}
    		if(m_mforce.getBooleanValue()){
    			command.add("-f");
    		}else{
    			if(new File (outfile).exists()){
    				throw new InvalidSettingsException("Outfile already exists. Please move/delete the file or enable outfile overwrite.");
    			}
    		}
    		if(m_msorted.getBooleanValue()){
    			command.add("-n");
    		}
    		if(m_mrgtag.getBooleanValue()){
    			command.add("-r");
    		}
    		if(m_muncompressed.getBooleanValue()){
    			command.add("-u");
    		}
    		if(m_usemhfile.getBooleanValue()){
    			command.add("-h "+IO.processFilePath(m_mhfile.getStringValue()));
    		}
    		if(m_usemregion.getBooleanValue()){
    			command.add("-R "+m_mregion.getStringValue());
    		}
    		
    		command.add(outfile);
    		
            Iterator<DataRow> it = inData[0].iterator();
            int size = 0;
            while(it.hasNext()){
            	command.add(it.next().getCell(colIndex).toString());
            	size++;
            } 
     		if(size<2){
     			throw new InvalidSettingsException("Insufficient number of input files. At least 2 bam files required for samtools cat.");
     		}
     		
    		OutCellType = "BAMCell";
    		
    	} 	else if (utility.equals("phase")) {

    		command.add("-b "+m_prefix.getStringValue());
    		command.add("-k "+m_blocklength.getIntValue());
    		command.add("-q "+m_hetphred.getIntValue());
    		command.add("-Q "+m_minqual.getIntValue());
    		command.add("-D "+m_maxdepth.getIntValue());
       		if(m_fixchimeras.getBooleanValue()){
       			command.add("-F");
       		}

       		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")+1) + "phase";
       		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
       		command.add(path2bamfile);
       		stdOut = true;
       		OutCellType = "BAMCell";
    	}
   	
    	
    	if(command.size()!=0) {
    		
    		if(!stdOut){	//No extra output file required
    			super.executeCommand(new String[]{StringUtils.join(command, " ")}, outfile, exec,new File(lockFile),null,outfile+".stdErr");
    			
    		}else{	//StdOut to outfile
    			super.executeCommand(new String[]{StringUtils.join(command, " ")}, outfile, exec, null, new File(lockFile), outfile,outfile+".stdErr", null, null, null);
    		}
    	}else{
    		throw new Exception("Something went wrong...Execution command is empty!");
    	}
   	
    	
    	DataColumnSpec dcs = null;
    	if(OutCellType.equals("BAMCell")){
    		dcs = new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec();
    	}else{
    		dcs = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					dcs}));
    	
    	FileCell[] c = new FileCell[]{ FileCellFactory.create(outfile)};
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
    	samtools_bin = IO.processFilePath(m_samtools.getStringValue());
    	ref_genome = IO.processFilePath(m_refseqfile.getStringValue());
   	
		if(CompatibilityChecker.inputFileNotOk(samtools_bin, false)) {
			throw new InvalidSettingsException("Set path to samtools binary!");
		}
       	
     	if(ref_genome.length() > 1) {
	     	if(!FileValidator.checkFastaFormat(ref_genome)){
	     		throw new InvalidSettingsException("Sequence file is not in fasta format!");
	     	}
     	}

    	DataColumnSpec dcs = null;
    	if(OutCellType.equals("BAMCell")){
    		dcs = new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec();
    	}else{
    		dcs = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();
    	}
    	    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					dcs})};
    }
}

