package de.helmholtz_muenchen.ibis.ngs.bcftools;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of Bcftools.
 * 
 *
 * @author Maximilian Hastreiter

 */
public class BcftoolsNodeModel extends HTExecutorNodeModel {
    /**
     * Input/output options:

       -A        keep all possible alternate alleles at variant sites
       -b        output BCF instead of VCF
       -D FILE   sequence dictionary for VCF->BCF conversion [null]
       -F        PL generated by r921 or before (which generate old ordering)
       -G        suppress all individual genotype information
       -l FILE   list of sites (chr pos) or regions (BED) to output [all sites]
       -L        calculate LD for adjacent sites
       -N        skip sites where REF is not A/C/G/T
       -Q        output the QCALL likelihood format
       -s FILE   list of samples to use [all samples]
       -S        input is VCF
       -u        uncompressed BCF output (force -b)

Consensus/variant calling options:

       -c        SNP calling (force -e)
       -d FLOAT  skip loci where less than FLOAT fraction of samples covered [0]
       -e        likelihood based analyses
       -g        call genotypes at variant sites (force -c)
       -i FLOAT  indel-to-substitution ratio [-1]
       -I        skip indels
       -p FLOAT  variant if P(ref|D)<FLOAT [0.5]
       -P STR    type of prior: full, cond2, flat [full]
       -t FLOAT  scaled substitution mutation rate [0.001]
       -T STR    constrained calling; STR can be: pair, trioauto, trioxd and trioxs (see manual) [null]
       -v        output potential variant sites only (force -c)

Contrast calling and association test options:

       -1 INT    number of group-1 samples [0]
       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [1]
       -U INT    number of permutations for association testing (effective with -1) [0]
       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [0.01]
     */
	
	/**
	 * Bcftools method chooser and infile
	 */
	public static final String CFGKEY_BCFMETHOD="bcfmethod";
	public static final String CFGKEY_INFILE="infile";
	public static final String CFGKEY_CATINFILE="catinfile";
	public static final String CFGKEY_PATH2BCFTOOLS="path2bcftools";
	public static final String CFGKEY_VCFSAMPLEHEADER="VCFSAMPLEHEADER";
	
	
	/**
	 * Bcftools call In/Out Options
	 */
	public static final String CFGKEY_KEEPALLELES="keepalleles";
	public static final String CFGKEY_OUTBCF="outbcf";
	public static final String CFGKEY_SEQDIC="seqdic";
	public static final String CFGKEY_PLGENERATE="plgenerate";
	public static final String CFGKEY_SURPESSGENOTYPE="surpressgenotype";
	public static final String CFGKEY_BEDFILE="bedfile";
	public static final String CFGKEY_CALCLD="calcld";
	public static final String CFGKEY_SKIPREF="skipref";
	public static final String CFGKEY_OUTQCALL="outqcall";
	public static final String CFGKEY_SAMPLELIST="samplelist";
	public static final String CFGKEY_OUTUNCOMPRESSEDBCF="outuncompressedbcf";

	//Checkboxes
	public static final String CFGKEY_IFBEDFILE="ifbedfile";
	public static final String CFGKEY_IFSEQDIC="ifseqdic";
	public static final String CFGKEY_IFSAMPLELIST="ifsamplelist";
	/**
	 * Bcftools view Consensus/variant calling options:
	 */
	public static final String CFGKEY_SNPCALLING="snpcalling";
	public static final String CFGKEY_SAMPLECOVERAGE="samplecoverage";
	public static final String CFGKEY_LIKELIHOODANA="likelihoodana";
	public static final String CFGKEY_CALLGENOTYPE="callgenotype";
	public static final String CFGKEY_INDELSUBRATIO="indelsubratio";
	public static final String CFGKEY_SKIPINDEL="skipindel";
	public static final String CFGKEY_VARIANTIF="variantif";
	public static final String CFGKEY_TYPEOFPRIOR="typeofprior";
	public static final String CFGKEY_MUTATIONRATE="mutationrate";
	public static final String CFGKEY_CONSTRAINEDCALLING="constrainedcalling";
	public static final String CFGKEY_OUTVARIANTSONLY="outvariantsonly";
	/**
	 * Contrast calling and association test options
	 */
	public static final String CFGKEY_NUMBEROFGRPSAMPLES="numberofgrpsamples";
	public static final String CFGKEY_POSTERIORICON="posterioricons";
	public static final String CFGKEY_NUMBEROFPERMUTAIONS="numberofpermuations";
	public static final String CFGKEY_ONLYPERMUTAIONS="onlypermuations";
	
	/**
	 * Default values for Bcftools view Consensus/variant calling options:
	 */
	public static final int DEFAULT_SAMPLECOVERAGE=0;
	public static final double DEFAULT_INDELSUBRATIO=-1.0;
	public static final double DEFAULT_VARIANTIF=0.5;
	public static final String DEFAULT_TYPEOFPRIOR="full";
	public static final double DEFAULT_MUTATIONRATE=0.001;
	/**
	 * Default values for Contrast calling and association test options
	 */
	public static final int DEFAULT_NUMBEROFGRPSAMPLES=0;
	public static final double DEFAULT_POSTERIORICON=1;
	public static final int DEFAULT_NUMBEROFPERMUTAIONS=0;
	public static final double DEFAULT_ONLYPERMUTAIONS=0.01;
	
	
	/**
	 * Bcf methods and inputfile
	 */
	private final SettingsModelString m_path2bcftools = new SettingsModelString(
			CFGKEY_PATH2BCFTOOLS, "");
	private final SettingsModelString m_bcfmethod = new SettingsModelString(
			CFGKEY_BCFMETHOD, "");
	private final SettingsModelString m_infile = new SettingsModelString(
			CFGKEY_INFILE, "");
	private final SettingsModelString m_catinfile = new SettingsModelString(
			CFGKEY_CATINFILE, "");
	private final SettingsModelString m_vcfsampleheader = new SettingsModelString(
			CFGKEY_VCFSAMPLEHEADER, "");
	
	
	/**
	 * Bcftools view In/Out Options models
	 **/
	private final SettingsModelBoolean m_keepallel = new SettingsModelBoolean(
			CFGKEY_KEEPALLELES, false);
	private final SettingsModelBoolean m_outbcf = new SettingsModelBoolean(
			CFGKEY_OUTBCF, false);
	private final SettingsModelString m_seqdic = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_SEQDIC,"");
	private final SettingsModelBoolean m_plgenerate = new SettingsModelBoolean(
			CFGKEY_PLGENERATE, false);
	private final SettingsModelBoolean m_surpressgenotype = new SettingsModelBoolean(
			CFGKEY_SURPESSGENOTYPE, false);
	private final SettingsModelString m_bedfile = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_BEDFILE,"");
	private final SettingsModelBoolean m_calcld = new SettingsModelBoolean(
			CFGKEY_CALCLD, false);
	private final SettingsModelBoolean m_skipref = new SettingsModelBoolean(
			CFGKEY_SKIPREF, false);
	private final SettingsModelBoolean m_outqcall = new SettingsModelBoolean(
			CFGKEY_OUTQCALL, false);
	private final SettingsModelString m_samplelist = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_SAMPLELIST,"");
	private final SettingsModelBoolean m_outuncompressedbcf = new SettingsModelBoolean(
			CFGKEY_OUTUNCOMPRESSEDBCF, false);
	//Checkboxes
	private final SettingsModelBoolean m_ifbedfile = new SettingsModelBoolean(
			CFGKEY_IFBEDFILE, false);
	private final SettingsModelBoolean m_ifseqdic = new SettingsModelBoolean(
			CFGKEY_IFSEQDIC, false);
	private final SettingsModelBoolean m_ifsamplelist = new SettingsModelBoolean(
			CFGKEY_IFSAMPLELIST, false);
	
	
	/**
	 * Bcftools view Consensus/variant calling options models
	 */
	private final SettingsModelBoolean m_snpcalling = new SettingsModelBoolean(
			CFGKEY_SNPCALLING, false);
	private final SettingsModelDoubleBounded m_samplecoverage = new SettingsModelDoubleBounded(
			CFGKEY_SAMPLECOVERAGE,DEFAULT_SAMPLECOVERAGE,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_lielihoodana = new SettingsModelBoolean(
			CFGKEY_LIKELIHOODANA, false);
	private final SettingsModelBoolean m_callgenotype = new SettingsModelBoolean(
			CFGKEY_CALLGENOTYPE, false);
	private final SettingsModelDoubleBounded m_indelsubratio = new SettingsModelDoubleBounded(
			CFGKEY_INDELSUBRATIO,DEFAULT_INDELSUBRATIO,-Double.MAX_VALUE,Double.MAX_VALUE);
	private final SettingsModelBoolean m_skipindel = new SettingsModelBoolean(
			CFGKEY_SKIPINDEL, false);
	private final SettingsModelDoubleBounded m_variantif = new SettingsModelDoubleBounded(
			CFGKEY_VARIANTIF,DEFAULT_VARIANTIF,0,Double.MAX_VALUE);
	private final SettingsModelString m_typeofprior = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_TYPEOFPRIOR,BcftoolsNodeModel.DEFAULT_TYPEOFPRIOR);
	private final SettingsModelDoubleBounded m_mutationrate = new SettingsModelDoubleBounded(
			CFGKEY_MUTATIONRATE,DEFAULT_MUTATIONRATE,0,Double.MAX_VALUE);
	private final SettingsModelString m_constrainedcalling = new SettingsModelString(
			BcftoolsNodeModel.CFGKEY_CONSTRAINEDCALLING,"");
	private final SettingsModelBoolean m_outvariantsonly = new SettingsModelBoolean(
			CFGKEY_OUTVARIANTSONLY, false);
	/**
	 * Contrast calling and association test models
	 */
	private final SettingsModelIntegerBounded m_numberofgrpsamples = new SettingsModelIntegerBounded(
			CFGKEY_NUMBEROFGRPSAMPLES,DEFAULT_NUMBEROFGRPSAMPLES,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_posterioricon = new SettingsModelDoubleBounded(
			CFGKEY_POSTERIORICON,DEFAULT_POSTERIORICON,0,Double.MAX_VALUE);
	private final SettingsModelIntegerBounded m_numberofpermuations = new SettingsModelIntegerBounded(
			CFGKEY_NUMBEROFPERMUTAIONS,DEFAULT_NUMBEROFPERMUTAIONS,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_onlypermutations = new SettingsModelDoubleBounded(
			CFGKEY_ONLYPERMUTAIONS,DEFAULT_ONLYPERMUTAIONS,0,Double.MAX_VALUE);
	
	
	private boolean OptionalPort=false;
	
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(BcftoolsNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2Bcftools";
	public static final String OUT_COL2 = "Path2BcftoolsOutput";
	
    /**
     * Constructor for the node model.
     */
    protected BcftoolsNodeModel() {
    
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
        
        m_bedfile.setEnabled(false);
        m_seqdic.setEnabled(false);
        m_samplelist.setEnabled(false);
        m_catinfile.setEnabled(false);
    	m_samplecoverage.setEnabled(false);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String infile="";
    	String path2bcftools="";
    	if(OptionalPort){
        	infile = inData[0].iterator().next().getCell(1).toString();
        	m_infile.setStringValue(infile);
        	path2bcftools=inData[0].iterator().next().getCell(0).toString();
        	path2bcftools=path2bcftools.substring(0, path2bcftools.lastIndexOf("/"));
        	path2bcftools+="/bcftools";
        	m_path2bcftools.setStringValue(path2bcftools);
    	}else{
        	infile =m_infile.getStringValue();
        	path2bcftools=m_path2bcftools.getStringValue();
    	}
      	 	
    	ArrayList<String> command = new ArrayList<String>();

    	String method = m_bcfmethod.getStringValue();
    	String outfile=m_infile.getStringValue();
    	command.add(path2bcftools+" "+method);
    	/**
    	 * About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.
         This command replaces the former "bcftools view" caller. Some of the original
         functionality has been temporarily lost in the process of transition to htslib,
         but will be added back on popular demand. The original calling model can be
         invoked with the -c option.
Usage:   bcftools call [options] <in.vcf.gz>

File format options:
   -o, --output <file>             write output to a file [standard output]
   -O, --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
   -r, --regions <region>          restrict to comma-separated list of regions
   -R, --regions-file <file>       restrict to regions listed in a file
   -s, --samples <list>            list of samples to include [all samples]
   -S, --samples-file <file>       PED file or a file with optional second column for ploidy (0, 1 or 2) [all samples]
   -t, --targets <region>          similar to -r but streams rather than index-jumps
   -T, --targets-file <file>       similar to -R but streams rather than index-jumps


    	 */
    	if(m_bcfmethod.getStringValue().equals("view")){
    	//Bcftools View !
    	/**
    	 *    Input/output options
    	 *    Input/output options:
   -A, --keep-alts                 keep all possible alternate alleles at variant sites
   -f, --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
   -g, --gvcf <minDP>              output gVCF blocks of homozygous REF calls
   -i, --insert-missed             output also sites missed by mpileup but present in -T
   -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)
   -V, --skip-variants <type>      skip indels/snps
   -v, --variants-only             output variant sites only

    	 */
    	if(m_keepallel.getBooleanValue()){
    		command.add("-A");
    	}
    
//    	if(m_seqdic.isEnabled()){command.add("-D "+m_seqdic.getStringValue());}
//    	if(m_plgenerate.getBooleanValue()){command.add(" -F");}
//    	if(m_surpressgenotype.getBooleanValue()){command.add(" -G");}
    	if(m_bedfile.isEnabled()){
    		command.add("-R "+m_bedfile.getStringValue());
    	}
//    	if(m_calcld.getBooleanValue()){command.add("-L");}
//    	if(m_skipref.getBooleanValue()){command.add("-N");}
//    	if(m_outqcall.getBooleanValue()){command.add("-Q");}
    	if(m_samplelist.isEnabled()){
    		command.add("-s "+m_samplelist.getStringValue());
    	}
    	
    	String fileformat = m_infile.getStringValue().substring(m_infile.getStringValue().lastIndexOf(".")+1);
    	
    	if(fileformat.equals("vcf")){
    		command.add("-O v");
    	}
    	if(m_outbcf.getBooleanValue()){
    		command.add("-O b");
    	}
    	if(m_outuncompressedbcf.getBooleanValue()){
    		command.add("-O u");
    	}
    	
    	/**
    	 *  Bcftools view Consensus/variant calling options models
    	 *  Consensus/variant calling options:
   -c, --consensus-caller          the original calling method (conflicts with -m)
   -C, --constrain <str>           one of: alleles, trio (see manual)
   -m, --multiallelic-caller       alternative model for multiallelic and rare-variant calling (conflicts with -c)
   -n, --novel-rate <float>,[...]  likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]
   -p, --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
   -P, --prior <float>             mutation rate (use bigger for greater sensitivity) [1.1e-3]
   -X, --chromosome-X              haploid output for male samples (requires PED file with -s)
   -Y, --chromosome-Y              haploid output for males and skips females (requires PED file with -s)

    	 */
    	if(m_snpcalling.getBooleanValue()){
    		command.add("-c");
    	}
//    	if(m_samplecoverage.isEnabled()){command.add("-d "+m_samplecoverage.getDoubleValue());}
//    	if(m_lielihoodana.getBooleanValue()){command.add("-e");}
//    	if(m_callgenotype.getBooleanValue()){command.add("-g");}
//    	if(m_indelsubratio.isEnabled()){command.add("-i "+m_indelsubratio.getDoubleValue());}
//    	if(m_skipindel.getBooleanValue()){command.add("-I");}
//    	command.add("-p "+m_variantif.getDoubleValue());
//    	command.add("-P "+m_typeofprior.getStringValue());
//    	command.add("-t "+m_mutationrate.getDoubleValue());
    	if(!m_constrainedcalling.getStringValue().equals("No constrains")){
    		command.add("-C "+m_constrainedcalling.getStringValue());
    	}
    	if(m_outvariantsonly.getBooleanValue()){
    		command.add("-v");
    	}
    	/**
    	 * Contrast calling and association test models
    	 */
//    	command.add("-1 "+m_numberofgrpsamples.getIntValue()); // this tag doesn't exist as of 28.4.15
//    	command.add("-C "+m_posterioricon.getDoubleValue());
//    	command.add("-U "+m_numberofpermuations.getIntValue());
//    	command.add("-X "+m_onlypermutations.getDoubleValue());
    	
    	/**
    	 * Infile & Outfile
    	 */
    	command.add(m_infile.getStringValue());

    	if(m_outbcf.getBooleanValue() || m_outuncompressedbcf.getBooleanValue()){
    		outfile+="_call.bcf";
    	}else{
    		outfile+="_call.vcf";
    	}
    	command.add("-o " +outfile);

    	}
    	else if(m_bcfmethod.getStringValue().equals("concat")){
            command.add(m_infile.getStringValue());
            command.add(m_catinfile.getStringValue());
        	String secfile =m_catinfile.getStringValue().substring(m_catinfile.getStringValue().lastIndexOf("/")+1,m_catinfile.getStringValue().length());
        	secfile=secfile.substring(0,secfile.lastIndexOf("."));
        	outfile+="_"+secfile+"_concat.bcf";
    	}
    	else if(m_bcfmethod.getStringValue().equals("index")){
            command.add(m_infile.getStringValue());
    	}
    	else if(m_bcfmethod.getStringValue().equals("reheader")){
            command.add("-s "+m_vcfsampleheader.getStringValue());
            command.add(m_infile.getStringValue());
            
            if(m_infile.getStringValue().endsWith(".vcf")){
            	outfile+="_reheader.vcf";
            }else if(m_infile.getStringValue().endsWith(".bcf")){
            	outfile+="_reheader.bcf";
            }else{
            	
            }
    	}

    	/**
    	 * Execute
    	 */
    	
    	String lockFile = outfile +"_"+m_bcfmethod.getStringValue()+SuccessfulRunChecker.LOCK_ENDING;
    	
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,null,LOGGER,outfile,outfile+".stdErr");
		super.executeCommand(new String[]{StringUtils.join(command, " ")},exec,new File(lockFile),outfile,outfile+".stdErr");
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2bcftools.substring(0,path2bcftools.lastIndexOf("/")+1)),
    			(FileCell) FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	pushFlowVariableString("BCFTOOLSOUT", outfile);
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
    		
    			//Check OptionalInputPort
    		try{
    			inSpecs[0].getNumColumns();
    			//Input table available
    			OptionalPort=true;
    			String[] colnames = inSpecs[0].getColumnNames();
    			if(colnames[0].equals("Path2SamTools")){//Path2Samtools available
    				m_path2bcftools.setEnabled(false);
    			}else{
//    				throw new InvalidSettingsException("First column name should be 'Path2SamTools' but it is" +inSpecs[0].getColumnNames()[0]);
    			}
//    			if(colnames[1].equals("Path2MpileupOutfile")){//Path2MpileupOutfile available
//    				m_infile.setEnabled(false);
//    			}else{
//    				throw new InvalidSettingsException("Column name of second column should be 'Path2MpileupOutfile' but it is" +inSpecs[0].getColumnNames()[1]);
//    			}
            	String fileformat = m_infile.getStringValue().substring(m_infile.getStringValue().lastIndexOf(".")+1);
    			if(fileformat.equals("vcf")&&m_outbcf.getBooleanValue()&&m_bcfmethod.getStringValue().equals("call")&&!(m_ifseqdic.getBooleanValue())){
    				throw new InvalidSettingsException("Input is VCF. Output is BCF. Please specify the sequence dictionary.");
    			}
    			if(!m_ifsamplelist.getBooleanValue()&&(m_constrainedcalling.getStringValue().equals("trioauto")||m_constrainedcalling.getStringValue().equals("trioxd")||m_constrainedcalling.getStringValue().equals("trioxs"))){
    				setWarningMessage("WARNING: You are using trio calling without a specified samplelist file. This may cause errors !");
    			}

    			return new DataTableSpec[]{new DataTableSpec(
    	    			new DataColumnSpec[]{
    	    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    	    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
    			
        	}catch(NullPointerException npe){
        			OptionalPort=false;
    				m_path2bcftools.setEnabled(true);
    				m_infile.setEnabled(true);
    	        	String fileformat = m_infile.getStringValue().substring(m_infile.getStringValue().lastIndexOf(".")+1);
    				if(fileformat.equals("vcf")&&m_outbcf.getBooleanValue()&&m_bcfmethod.getStringValue().equals("call")&&!(m_ifseqdic.getBooleanValue())){
    					throw new InvalidSettingsException("Input is VCF. Output is BCF. Please specify the sequence dictionary.");
    				}
    				if(!m_ifsamplelist.getBooleanValue()&&(m_constrainedcalling.getStringValue().equals("trioauto")||m_constrainedcalling.getStringValue().equals("trioxd")||m_constrainedcalling.getStringValue().equals("trioxs"))){
    					setWarningMessage("WARNING: You are using trio calling without a specified samplelist file. This may cause errors !");
    				}
    		    	
    		        return new DataTableSpec[]{new DataTableSpec(
    		    			new DataColumnSpec[]{
    		    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    		    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
        	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	 m_bcfmethod.saveSettingsTo(settings);
    	 m_infile.saveSettingsTo(settings);
    	 m_catinfile.saveSettingsTo(settings);
    	 m_path2bcftools.saveSettingsTo(settings);
    	 m_vcfsampleheader.saveSettingsTo(settings);
    	
    	 m_bedfile.saveSettingsTo(settings);
         m_calcld.saveSettingsTo(settings);
         m_callgenotype.saveSettingsTo(settings);
         m_constrainedcalling.saveSettingsTo(settings);
         m_indelsubratio.saveSettingsTo(settings);
         m_keepallel.saveSettingsTo(settings);
         m_lielihoodana.saveSettingsTo(settings);
         m_mutationrate.saveSettingsTo(settings);
         m_numberofgrpsamples.saveSettingsTo(settings);
         m_numberofpermuations.saveSettingsTo(settings);
         m_onlypermutations.saveSettingsTo(settings);
         m_outbcf.saveSettingsTo(settings);
         m_outqcall.saveSettingsTo(settings);
         m_outuncompressedbcf.saveSettingsTo(settings);
         m_outvariantsonly.saveSettingsTo(settings);
         m_plgenerate.saveSettingsTo(settings);
         m_posterioricon.saveSettingsTo(settings);
         m_samplecoverage.saveSettingsTo(settings);
         m_samplelist.saveSettingsTo(settings);
         m_seqdic.saveSettingsTo(settings);
         m_skipindel.saveSettingsTo(settings);
         m_skipref.saveSettingsTo(settings);
         m_snpcalling.saveSettingsTo(settings);
         m_surpressgenotype.saveSettingsTo(settings);
         m_typeofprior.saveSettingsTo(settings);
         m_variantif.saveSettingsTo(settings);
         
         m_ifbedfile.saveSettingsTo(settings);
         m_ifsamplelist.saveSettingsTo(settings);
         m_ifseqdic.saveSettingsTo(settings);
         
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_bcfmethod.loadSettingsFrom(settings);
    	m_infile.loadSettingsFrom(settings);
    	m_catinfile.loadSettingsFrom(settings);
   	 	m_path2bcftools.loadSettingsFrom(settings);
   	 	m_vcfsampleheader.loadSettingsFrom(settings);
    	
        m_bedfile.loadSettingsFrom(settings);
        m_calcld.loadSettingsFrom(settings);
        m_callgenotype.loadSettingsFrom(settings);
        m_constrainedcalling.loadSettingsFrom(settings);
        m_indelsubratio.loadSettingsFrom(settings);
        m_keepallel.loadSettingsFrom(settings);
        m_lielihoodana.loadSettingsFrom(settings);
        m_mutationrate.loadSettingsFrom(settings);
        m_numberofgrpsamples.loadSettingsFrom(settings);
        m_numberofpermuations.loadSettingsFrom(settings);
        m_onlypermutations.loadSettingsFrom(settings);
        m_outbcf.loadSettingsFrom(settings);
        m_outqcall.loadSettingsFrom(settings);
        m_outuncompressedbcf.loadSettingsFrom(settings);
        m_outvariantsonly.loadSettingsFrom(settings);
        m_plgenerate.loadSettingsFrom(settings);
        m_posterioricon.loadSettingsFrom(settings);
        m_samplecoverage.loadSettingsFrom(settings);
        m_samplelist.loadSettingsFrom(settings);
        m_seqdic.loadSettingsFrom(settings);
        m_skipindel.loadSettingsFrom(settings);
        m_skipref.loadSettingsFrom(settings);
        m_snpcalling.loadSettingsFrom(settings);
        m_surpressgenotype.loadSettingsFrom(settings);
        m_typeofprior.loadSettingsFrom(settings);
        m_variantif.loadSettingsFrom(settings);
        m_ifbedfile.loadSettingsFrom(settings);
        m_ifsamplelist.loadSettingsFrom(settings);
        m_ifseqdic.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
   	 	m_bcfmethod.validateSettings(settings);
   	 	m_infile.validateSettings(settings);
   	 	m_catinfile.validateSettings(settings);
   	 	m_path2bcftools.validateSettings(settings);
   	 	m_vcfsampleheader.validateSettings(settings);
    	
    	m_bedfile.validateSettings(settings);
        m_calcld.validateSettings(settings);
        m_callgenotype.validateSettings(settings);
        m_constrainedcalling.validateSettings(settings);
        m_indelsubratio.validateSettings(settings);
        m_keepallel.validateSettings(settings);
        m_lielihoodana.validateSettings(settings);
        m_mutationrate.validateSettings(settings);
        m_numberofgrpsamples.validateSettings(settings);
        m_numberofpermuations.validateSettings(settings);
        m_onlypermutations.validateSettings(settings);
        m_outbcf.validateSettings(settings);
        m_outqcall.validateSettings(settings);
        m_outuncompressedbcf.validateSettings(settings);
        m_outvariantsonly.validateSettings(settings);
        m_plgenerate.validateSettings(settings);
        m_posterioricon.validateSettings(settings);
        m_samplecoverage.validateSettings(settings);
        m_samplelist.validateSettings(settings);
        m_seqdic.validateSettings(settings);
        m_skipindel.validateSettings(settings);
        m_skipref.validateSettings(settings);
        m_snpcalling.validateSettings(settings);
        m_surpressgenotype.validateSettings(settings);
        m_typeofprior.validateSettings(settings);
        m_variantif.validateSettings(settings);
        m_ifbedfile.validateSettings(settings);
        m_ifsamplelist.validateSettings(settings);
        m_ifseqdic.validateSettings(settings);
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

