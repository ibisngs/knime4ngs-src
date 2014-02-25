package de.helmholtz_muenchen.ibis.ngs.bamsamconverter;

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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of BAMSAMConverter.
 * 
 *
 * @author Maximilian Hastreiter
 * @author Sebastian Kopetzky
 */
public class BAMSAMConverterNodeModel extends NodeModel {
    
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BAMSAMConverterNodeModel.class);
	
	public static final String CFGKEY_PATH2SAMTOOLS="path2samtools";
	public static final String CFGKEY_METHOD="method";
	public static final String CFGKEY_INFILE="infile";
	
	public static final String CFGKEY_printsamheader="printsamheader";
	public static final String CFGKEY_printsamheaderonly="printsamheaderonly";
	public static final String CFGKEY_uncompressedbam="uncompressedbam";
	public static final String CFGKEY_fastcompression="fastcompression";
	public static final String CFGKEY_outhexflag="outhexflag";
	public static final String CFGKEY_outstringflag="outstringflag";
	public static final String CFGKEY_printmatchingrecordsonly="printmatchingrecordsonly";
	public static final String CFGKEY_bedfile="bedfile";
	public static final String CFGKEY_refnamelist="refnamelist";
	public static final String CFGKEY_refseqfile="refseqfile";
	public static final String CFGKEY_listofreads="listofreads";
	public static final String CFGKEY_requiredflag="requiredflag";
	public static final String CFGKEY_filteringflag="filteringflag";
	public static final String CFGKEY_minmapqual="minmapqual";
	public static final String CFGKEY_onlylibraryreads="onlylibraryreads";
	public static final String CFGKEY_onlygroupreads="onlygroupreads";
	
	//Checkboxes

	public static final String CFGKEY_USEBEDFILE="usebedfile";
	public static final String CFGKEY_USEREFNAMELIST="userefnamelist";
	public static final String CFGKEY_USEREFSEQFILE="userefseqfile";
	public static final String CFGKEY_USELISTOFREADS="uselistofreads";
	public static final String CFGKEY_USEREQUIREDFLAG="userequiredflag";
	public static final String CFGKEY_USEFILTERINGFLAG="usefilteringflag";
	public static final String CFGKEY_USEONLYLIBRARYREADS="useonlylibraryreads";
	public static final String CFGKEY_USEONLYGROUPREADS="useonlygroupreads";
	public static final String CFGKEY_CHECKBAMINDEX="checkbamindex";
	
	private final SettingsModelString m_path2samtools = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_PATH2SAMTOOLS,"");
	private final SettingsModelString m_method = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_METHOD,"");
	private final SettingsModelString m_infile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_INFILE,"");
	private final SettingsModelBoolean m_printsamheader = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printsamheader, true);
	private final SettingsModelBoolean m_printsamheaderonly = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printsamheaderonly, false);
	private final SettingsModelBoolean m_uncompressedbam = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_uncompressedbam, false);
	private final SettingsModelBoolean m_fastcompression = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_fastcompression, false);
	private final SettingsModelBoolean m_outhexflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_outhexflag, false);
	private final SettingsModelBoolean m_outstringflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_outstringflag, false);
	private final SettingsModelBoolean m_printmatchingrecordsonly = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printmatchingrecordsonly, false);
	private final SettingsModelString m_bedfile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_bedfile, "");
	private final SettingsModelString m_refnamelist = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_refnamelist, "");
	private final SettingsModelString m_refseqfile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_refseqfile, "");
	private final SettingsModelString m_listofreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_listofreads, "");	
	private final SettingsModelIntegerBounded m_requiredflag = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_requiredflag,0,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_filteringflag = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_filteringflag,0,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minmapqual = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_minmapqual,0,0,Integer.MAX_VALUE);
	private final SettingsModelString m_onlylibraryreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_onlylibraryreads, "");
	private final SettingsModelString m_onlygroupreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_onlygroupreads, "");
	
	//Checkboxes
	private final SettingsModelBoolean m_usebedfile = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEBEDFILE, false);
	private final SettingsModelBoolean m_userefnamelist = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREFNAMELIST, false);
	private final SettingsModelBoolean m_userefseqfile = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREFSEQFILE, false);
	private final SettingsModelBoolean m_userequiredflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREQUIREDFLAG, false);
	private final SettingsModelBoolean m_usefilteringflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEFILTERINGFLAG, false);
	private final SettingsModelBoolean m_useonlylibraryreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEONLYLIBRARYREADS, false);
	private final SettingsModelBoolean m_useonlygroupreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEONLYGROUPREADS, false);
	private final SettingsModelBoolean m_uselistofreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USELISTOFREADS, false);
	private final SettingsModelBoolean m_checkbamindex = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_CHECKBAMINDEX, true);
	
	private boolean OptionalPort=false;
	private boolean INSAM0=false;
	private boolean INSAM1=false;
	private boolean INBAM=false;
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SamTools";
	public static String OUT_COL2 = "";
	
    /**
     * Constructor for the node model.
     */
    protected BAMSAMConverterNodeModel() {
    
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
        
        m_bedfile.setEnabled(false);
        m_refnamelist.setEnabled(false);
        m_refseqfile.setEnabled(false);
        m_listofreads.setEnabled(false);
        m_requiredflag.setEnabled(false);
        m_filteringflag.setEnabled(false);
        m_onlygroupreads.setEnabled(false);
        m_onlylibraryreads.setEnabled(false);
        m_printsamheader.setEnabled(false);
        m_printsamheaderonly.setEnabled(false);        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String infile="";
    	String path2samtools="";
    	if(OptionalPort){
    		if(INSAM0){
            	infile = inData[0].iterator().next().getCell(0).toString();
        		path2samtools=m_path2samtools.getStringValue();
        		m_infile.setStringValue(infile);
            	pushFlowVariableString("Path2seqFile", inData[0].iterator().next().getCell(1).toString());
    		}else if(INSAM1){
            	infile = inData[0].iterator().next().getCell(1).toString();
            	path2samtools=inData[0].iterator().next().getCell(0).toString();
        		m_infile.setStringValue(infile);
        		m_path2samtools.setStringValue(path2samtools);
            	pushFlowVariableString("Path2seqFile", inData[0].iterator().next().getCell(1).toString());
    		}else if(INBAM){
            	infile = inData[0].iterator().next().getCell(1).toString();	
            	path2samtools=inData[0].iterator().next().getCell(0).toString();
            	m_infile.setStringValue(infile);
            	m_path2samtools.setStringValue(path2samtools);
            	pushFlowVariableString("Path2seqFile", inData[0].iterator().next().getCell(2).toString());
    		}
    	}else{
    		infile=m_infile.getStringValue();
    		path2samtools=m_path2samtools.getStringValue();
        	pushFlowVariableString("Path2seqFile","No reference file");
    	}
    	
    	
    	/**Initialize logfile**/
    	String folder = infile;
    	folder = folder.substring(0,folder.lastIndexOf("/")+1);
    	String logfile = folder +"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("BAMSAMConverter"));
    	/**logfile initialized**/
    	
    	String outputfile="";
    	if(m_method.getStringValue().equals("SAM->BAM")){
        outputfile=infile.substring(0,infile.lastIndexOf("."))+".bam";
    	}else if(m_method.getStringValue().equals("BAM->SAM")){
    	outputfile=infile.substring(0,infile.lastIndexOf("."))+".sam";
    	}

    	
    	ArrayList<String> command = new ArrayList<String>();
    	command.add(path2samtools+" view");
    	if(m_method.getStringValue().equals("SAM->BAM")){command.add("-b");}
    	if(m_printsamheader.getBooleanValue()&&m_printsamheader.isEnabled()){command.add("-h");}
    	if(m_printsamheaderonly.getBooleanValue()&&m_printsamheaderonly.isEnabled()){command.add("-H");}
    	if(m_method.getStringValue().equals("SAM->BAM")){command.add("-S");}
    	if(m_uncompressedbam.getBooleanValue()&&m_uncompressedbam.isEnabled()){command.add("-u");}
    	if(m_fastcompression.getBooleanValue()&&m_fastcompression.isEnabled()){command.add("-1");}
    	if(m_outhexflag.getBooleanValue()&&m_outhexflag.isEnabled()){command.add("-x");}
    	if(m_outstringflag.getBooleanValue()&&m_outstringflag.isEnabled()){command.add("-X");}
    	if(m_printmatchingrecordsonly.getBooleanValue()&&m_printmatchingrecordsonly.isEnabled()){command.add("-c");}
    	if(m_usebedfile.getBooleanValue()&&m_usebedfile.isEnabled()){command.add("-L "+m_bedfile.getStringValue());}
    	if(m_userefnamelist.getBooleanValue()&&m_userefnamelist.isEnabled()){command.add("-t "+m_refnamelist.getStringValue());}
    	if(m_userefseqfile.getBooleanValue()&&m_userefseqfile.isEnabled()){command.add("-T "+m_refseqfile.getStringValue());}
    	command.add("-o "+outputfile);
    	if(m_uselistofreads.getBooleanValue()&&m_uselistofreads.isEnabled()){command.add("-R "+m_listofreads.getStringValue());}
    	if(m_userequiredflag.getBooleanValue()&&m_userequiredflag.isEnabled()){command.add("-f "+m_requiredflag.getIntValue());}
    	if(m_usefilteringflag.getBooleanValue()&&m_usefilteringflag.isEnabled()){command.add("-F "+m_filteringflag.getIntValue());}
    	command.add("-q "+m_minmapqual.getIntValue());
    	if(m_useonlylibraryreads.getBooleanValue()&&m_useonlylibraryreads.isEnabled()){command.add("-l "+m_onlylibraryreads.getStringValue());}
    	if(m_useonlygroupreads.getBooleanValue()&&m_useonlygroupreads.isEnabled()){command.add("-s "+m_onlygroupreads.getStringValue());}
    	command.add(infile);

		Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
	    logBuffer.append(ShowOutput.getNodeEndTime());
	    ShowOutput.writeLogFile(logBuffer);
		
    	
	    /**
	     * Converting finished, now sort and index if required
	     */
	    
	    
    	if(m_method.getStringValue().equals("SAM->BAM")){
    		command = new ArrayList<String>();	//clear list
    		
		    	//Sort BAM file
		        LOGGER.info("BAM Sort Process");
		        String path2bamfilesorted = outputfile.substring(0,outputfile.length()-4)+"_sorted";
		        command.add(path2samtools+" sort");
		        command.add(outputfile);
		        command.add(path2bamfilesorted);
	
		        Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
				
		      	
	        if(m_checkbamindex.getBooleanValue()){
	            //Index BAM file 
	        	command = new ArrayList<String>();	//clear list
	        	path2bamfilesorted += ".bam";
	        	LOGGER.info("BAM Index Process");
	        	command.add(path2samtools+" index");
	        	command.add(path2bamfilesorted);
	        	
		        Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
	    		
	
	        }else{
	        	LOGGER.info("Skipping BAM Index Process");
	        }
	        
	    	outputfile=path2bamfilesorted;
    	}
    	
    	/**
    	 * OUTPUT
    	 */
    	
    	if(m_method.getStringValue().equals("SAM->BAM")){
    		OUT_COL2 = "Path2BAMFile";
    	}else{
			OUT_COL2 = "Path2SAMFile";
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2samtools),
    			(FileCell) FileCellFactory.create(outputfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	pushFlowVariableString("MPILEUPINFILE", outputfile);

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
    	try{
        	inSpecs[0].getNumColumns();
        	//In-Port available
        	OptionalPort=true;
        	
        	/**
        	 * Check if BAMSAMINFILE variable exists.
        	 */
        	String infile="";
        	if(getAvailableInputFlowVariables().containsKey("BAMSAMINFILE")){
        		infile=getAvailableInputFlowVariables().get("BAMSAMINFILE").getStringValue();
        		}else{
        			throw new InvalidSettingsException("Infile variable not available. Wrong input node ! Requires Bfast,BWA,Bowtie2,Segemehl,BAMLoader or SAMLoader.");
        		}
        	String[] columnnames=inSpecs[0].getColumnNames();
        	if(columnnames[0].equals("Path2SamTools")){//BAMLoader || BAMSAMConverter
        		m_path2samtools.setEnabled(false);
        		if(columnnames[1].equals("Path2BAMFile")){
            		INBAM=true;            		
            		/**
            		 * Check if Outfile exists
            		 */		
            			String outfile=infile.substring(0,infile.lastIndexOf("."))+".sam";
                    	File outpath = new File(outfile);
                    	if(outpath.exists()){
                    		setWarningMessage(outfile+" already exists ! Please rename or move to other directory.");
                    	}
            		
                    //Disable fields for BAM->SAM                  	
        			m_infile.setEnabled(false);
        			m_method.setStringValue("BAM->SAM");
        			m_method.setEnabled(false);       			
        			m_uncompressedbam.setEnabled(false);
        			m_uncompressedbam.setBooleanValue(false);
        			m_fastcompression.setEnabled(false);
        			m_fastcompression.setBooleanValue(false);
        			m_checkbamindex.setEnabled(false);
        			m_checkbamindex.setBooleanValue(false);
        			m_userefnamelist.setEnabled(false);
        			m_userefseqfile.setEnabled(false);
        			m_printsamheader.setEnabled(true);
        			m_printsamheaderonly.setEnabled(true);      			
        			
        	    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SamTools", FileCell.TYPE);
        	    	DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2SAMFile", FileCell.TYPE);
        	    	DataTableSpec out = new DataTableSpec(new DataColumnSpec[]{col1.createSpec(),col2.createSpec()});
        	        return new DataTableSpec[]{out};		
        		}else if(columnnames[1].equals("Path2SAMFile")){
            		INSAM1=true;
            		/**
            		 * Check if Outfile exists
            		 */
                    	String outfile=infile.substring(0,infile.lastIndexOf("."))+".bam";
                    	File outpath = new File(outfile);
                    	if(outpath.exists()){
                    		setWarningMessage(outfile+" already exists ! Please rename or move to other directory.");
                    	}
          			//Disable/Enable fields for SAM->BAM
            		m_infile.setEnabled(false);
            		m_method.setStringValue("SAM->BAM");
            		m_method.setEnabled(false);
        			m_uncompressedbam.setEnabled(true);
        			m_fastcompression.setEnabled(true);
        			m_checkbamindex.setEnabled(true);
        			m_userefnamelist.setEnabled(true);
        			m_userefseqfile.setEnabled(true);
        			m_printsamheader.setEnabled(false);
        			m_printsamheaderonly.setEnabled(false);
            		
        	    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SamTools", FileCell.TYPE);
        	    	DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2BAMFile", FileCell.TYPE);
        	        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
        	    	DataTableSpec out = new DataTableSpec(cols);
        	        return new DataTableSpec[]{out};	
        		}else{
        			throw new InvalidSettingsException("Second column name should be 'Path2BAMFile or Path2SAMFile' but it is "+inSpecs[0].getColumnNames()[1]);
        		}
        	}else if(columnnames[0].equals("Path2SAMFile")){//SAMLoader || AlignerS
        		INSAM0=true;
        		/**
        		 * Check if Outfile exists
        		 */
                	String outfile=infile.substring(0,infile.lastIndexOf("."))+".bam";
                	File outpath = new File(outfile);
                	if(outpath.exists()){
                		setWarningMessage(outfile+" already exists ! Please rename or move to other directory.");
                	}
        		 	
                //Disable/Enable fields	for SAM->BAM
        		m_infile.setEnabled(false);
        		m_method.setStringValue("SAM->BAM");
        		m_method.setEnabled(false);
    			m_uncompressedbam.setEnabled(true);
    			m_fastcompression.setEnabled(true);
    			m_checkbamindex.setEnabled(true);
    			m_userefnamelist.setEnabled(true);
    			m_userefseqfile.setEnabled(true);
    			m_printsamheader.setEnabled(false);
    			m_printsamheaderonly.setEnabled(false);
        		
    	    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SamTools", FileCell.TYPE);
    	    	DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2BAMFile", FileCell.TYPE);
    	        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
    	    	DataTableSpec out = new DataTableSpec(cols);
    	    	
    	        return new DataTableSpec[]{out};	
        	}else{
        		throw new InvalidSettingsException("First column name should be 'Path2SamTools'(BAM->SAM) or 'SAM file'(SAM->BAM) but it is "+inSpecs[0].getColumnNames()[0]);
        	}

        	}catch(NullPointerException npe){
        		OptionalPort=false;
        		INSAM0=false;
        		INSAM1=false;
        		INBAM=false;
        		
    			//Disable/Enable fields for stand-alone
        		m_infile.setEnabled(true);
        		m_path2samtools.setEnabled(true);
        		m_method.setEnabled(true);
    			m_uncompressedbam.setEnabled(true);
    			m_fastcompression.setEnabled(true);
    			m_checkbamindex.setEnabled(true);
    			m_userefnamelist.setEnabled(true);
    			m_userefseqfile.setEnabled(true);
    			m_printsamheader.setEnabled(false);
    			m_printsamheaderonly.setEnabled(false);
    			
            	String conversion=m_method.getStringValue();
        		/**
        		 * Check if Outfile exists
        		 */
        			String outfile="";
                	if(conversion.equals("SAM->BAM")&&m_infile.getStringValue().length()>0){
                		outfile=m_infile.getStringValue().substring(0,m_infile.getStringValue().lastIndexOf("."))+".bam";
                	}else if(m_infile.getStringValue().length()>0&&conversion.equals("BAM->SAM")){
                		outfile=m_infile.getStringValue().substring(0,m_infile.getStringValue().lastIndexOf("."))+".sam";
                	} 			
                	File outpath = new File(outfile);
                	if(outpath.exists()){
                		setWarningMessage(outfile+" already exists ! Please rename or move to other directory.");
                	}
        		
        		
            	String file=m_infile.getStringValue();
            	String format=file.substring(file.lastIndexOf(".")+1);
            	if(conversion.equals("SAM->BAM")&&format.equals("bam")){
                    throw new InvalidSettingsException("SAM->BAM conversion requires SAM input file !");
            	}else if(conversion.equals("BAM->SAM")&&format.equals("sam")){
                    throw new InvalidSettingsException("BAM->SAM conversion requires BAM input file !");
            	}
    	    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SamTools", FileCell.TYPE);
    	    	DataColumnSpecCreator col2;
    	    	if(conversion.equals("SAM->BAM")){
    	    		col2 = new DataColumnSpecCreator("Path2BAMFile", FileCell.TYPE);
    	    	}else{
    	    		col2 = new DataColumnSpecCreator("Path2SAMFile", FileCell.TYPE);
    	    	}
    	    	DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
    	    	DataTableSpec out = new DataTableSpec(cols);
    	    	

    	        return new DataTableSpec[]{out};	
        	}    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_bedfile.saveSettingsTo(settings);
         m_fastcompression.saveSettingsTo(settings);
         m_filteringflag.saveSettingsTo(settings);
         m_infile.saveSettingsTo(settings);
         m_listofreads.saveSettingsTo(settings);
         m_minmapqual.saveSettingsTo(settings);
         m_onlygroupreads.saveSettingsTo(settings);
         m_onlylibraryreads.saveSettingsTo(settings);
         m_outhexflag.saveSettingsTo(settings);
         m_outstringflag.saveSettingsTo(settings);
         m_printmatchingrecordsonly.saveSettingsTo(settings);
         m_printsamheader.saveSettingsTo(settings);
         m_printsamheaderonly.saveSettingsTo(settings);
         m_refnamelist.saveSettingsTo(settings);
         m_refseqfile.saveSettingsTo(settings);
         m_requiredflag.saveSettingsTo(settings);
         m_uncompressedbam.saveSettingsTo(settings);
         m_usebedfile.saveSettingsTo(settings);
         m_usefilteringflag.saveSettingsTo(settings);
         m_useonlygroupreads.saveSettingsTo(settings);
         m_useonlylibraryreads.saveSettingsTo(settings);
         m_userefseqfile.saveSettingsTo(settings);
         m_userequiredflag.saveSettingsTo(settings);
         m_userefnamelist.saveSettingsTo(settings);
         m_path2samtools.saveSettingsTo(settings);
         m_method.saveSettingsTo(settings);
         m_uselistofreads.saveSettingsTo(settings);
         m_checkbamindex.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_bedfile.loadSettingsFrom(settings);
        m_fastcompression.loadSettingsFrom(settings);
        m_filteringflag.loadSettingsFrom(settings);
        m_infile.loadSettingsFrom(settings);
        m_listofreads.loadSettingsFrom(settings);
        m_minmapqual.loadSettingsFrom(settings);
        m_onlygroupreads.loadSettingsFrom(settings);
        m_onlylibraryreads.loadSettingsFrom(settings);
        m_outhexflag.loadSettingsFrom(settings);
        m_outstringflag.loadSettingsFrom(settings);
        m_printmatchingrecordsonly.loadSettingsFrom(settings);
        m_printsamheader.loadSettingsFrom(settings);
        m_printsamheaderonly.loadSettingsFrom(settings);
        m_refnamelist.loadSettingsFrom(settings);
        m_refseqfile.loadSettingsFrom(settings);
        m_requiredflag.loadSettingsFrom(settings);
        m_uncompressedbam.loadSettingsFrom(settings);
        m_usebedfile.loadSettingsFrom(settings);
        m_usefilteringflag.loadSettingsFrom(settings);
        m_useonlygroupreads.loadSettingsFrom(settings);
        m_useonlylibraryreads.loadSettingsFrom(settings);
        m_userefseqfile.loadSettingsFrom(settings);
        m_userequiredflag.loadSettingsFrom(settings);
        m_userefnamelist.loadSettingsFrom(settings);
        m_path2samtools.loadSettingsFrom(settings);
        m_method.loadSettingsFrom(settings);
        m_uselistofreads.loadSettingsFrom(settings);
        m_checkbamindex.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_bedfile.validateSettings(settings);
        m_fastcompression.validateSettings(settings);
        m_filteringflag.validateSettings(settings);
        m_infile.validateSettings(settings);
        m_listofreads.validateSettings(settings);
        m_minmapqual.validateSettings(settings);
        m_onlygroupreads.validateSettings(settings);
        m_onlylibraryreads.validateSettings(settings);
        m_outhexflag.validateSettings(settings);
        m_outstringflag.validateSettings(settings);
        m_printmatchingrecordsonly.validateSettings(settings);
        m_printsamheader.validateSettings(settings);
        m_printsamheaderonly.validateSettings(settings);
        m_refnamelist.validateSettings(settings);
        m_refseqfile.validateSettings(settings);
        m_requiredflag.validateSettings(settings);
        m_uncompressedbam.validateSettings(settings);
        m_usebedfile.validateSettings(settings);
        m_usefilteringflag.validateSettings(settings);
        m_useonlygroupreads.validateSettings(settings);
        m_useonlylibraryreads.validateSettings(settings);
        m_userefseqfile.validateSettings(settings);
        m_userequiredflag.validateSettings(settings);
        m_userefnamelist.validateSettings(settings);
        m_path2samtools.validateSettings(settings);
        m_method.validateSettings(settings);
        m_uselistofreads.validateSettings(settings);
        m_checkbamindex.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}

