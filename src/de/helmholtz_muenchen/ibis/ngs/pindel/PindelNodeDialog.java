package de.helmholtz_muenchen.ibis.ngs.pindel;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentDate;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentNumberEdit;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDate;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Pindel" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class PindelNodeDialog extends DefaultNodeSettingsPane {	
	
	/*
	 * executables
	 * String path to pindel
	 * 
	 * interval to analyse
	 * String interval default (ALL) (-c)
	 * (21 -> only chromosome 21, chromosome name has to match chromsome name in reference sequence
	 * 21:5,000,000-6,000,000 ->all bases in the interval on the chromosome)
	 * 
	 * pindel config file
	 * boolean create pindel config file (only previous node Picard insert size metrics)
	 * String path to pindel config file
	 * 
	 * int number of threads (-T)
	 * int window size : does analysis for x mio bases (Mb, megabases) -> influences RAM used and runtime
	 * 
	 * sensitivity/ selectivity
	 * int min num matched bases (30)
	 * int additional mismatch (1) -> only map part of read when no other position with the specified number of mismatches
	 * int min perfect match around BP (3) -> when mapping parts of reads
	 * float sequencing error rate (0.05)
	 * float maximum allowed mismatch rate for considering a reads (0.1)
	 * 
	 * boolean create vcf output -> pindel_SI, pindel_D 
	 * String path to pindel2vcf
	 *  
	 */
	
	// minimum coverage reads (10)
	// heterozygosity threshold (0.2)
	// homozygosity threshold (0.8)
	// minimum reads supporting the variant (1)
	// event supported on both strands? (boolean)
	// minimum size of variant (1)
	// maximum size of variant (infinite)
	// GATK compatible (boolean)
	
	final SettingsModelString pindel = new SettingsModelString(PindelNodeModel.CFGKEY_PINDEL, PindelNodeModel.DEF_PINDEL);;
	final SettingsModelBoolean interval = new SettingsModelBoolean(PindelNodeModel.CFGKEY_INTERVAL, PindelNodeModel.DEF_INTERVAL);;
	final SettingsModelString chrom= new SettingsModelString(PindelNodeModel.CFGKEY_CHROM, PindelNodeModel.DEF_CHROM);
	final SettingsModelIntegerBounded start = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_START, PindelNodeModel.DEF_START, PindelNodeModel.MIN_START, PindelNodeModel.MAX_START);
	final SettingsModelIntegerBounded end  = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_END, PindelNodeModel.DEF_END, PindelNodeModel.MIN_END, PindelNodeModel.MAX_END);
	final SettingsModelString config_file = new SettingsModelString(PindelNodeModel.CFGKEY_CONFIG_FILE, PindelNodeModel.DEF_CONFIG_FILE);
	final SettingsModelBoolean create_config  = new SettingsModelBoolean(PindelNodeModel.CFGKEY_CREATE_CONFIG, PindelNodeModel.DEF_CREATE_CONFIG);
	final SettingsModelBoolean vcf_out= new SettingsModelBoolean(PindelNodeModel.CFGKEY_VCF_OUT, PindelNodeModel.DEF_VCF_OUT);
	final SettingsModelString pindel2vcf = new SettingsModelString(PindelNodeModel.CFGKEY_VCF2PINDEL, PindelNodeModel.DEF_VCF2PINDEL);	
	final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_THREADS, PindelNodeModel.DEF_THREADS, PindelNodeModel.MIN_THREADS, PindelNodeModel.MAX_THREADS);
	final SettingsModelIntegerBounded bin_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_BIN_SIZE, PindelNodeModel.DEF_BIN_SIZE, PindelNodeModel.MIN_BIN_SIZE, PindelNodeModel.MAX_BIN_SIZE);
	
	// pindel params
	final SettingsModelInteger min_match_bases = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_MATCH_BASES, PindelNodeModel.DEF_MIN_MATCH_BASES, PindelNodeModel.MIN_MIN_MATCH_BASES, PindelNodeModel.MAX_MIN_MATCH_BASES);
	final SettingsModelInteger additional_mismatch = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_ADDITIONAL_MISMATCH, PindelNodeModel.DEF_ADDITIONAL_MISMATCH, PindelNodeModel.MIN_ADDITIONAL_MISMATCH, PindelNodeModel.MAX_ADDITIONAL_MISMATCH);
	final SettingsModelInteger min_match_breakpoint = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_MATCH_BP, PindelNodeModel.DEF_MIN_MATCH_BP, PindelNodeModel.MIN_MIN_MATCH_BP, PindelNodeModel.MAX_MIN_MATCH_BP);
	final SettingsModelDoubleBounded seq_err = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_SEQ_ERR, PindelNodeModel.DEF_SEQ_ERROR, PindelNodeModel.MIN_SEQ_ERROR, PindelNodeModel.MAX_SEQ_ERROR);
	final SettingsModelDoubleBounded max_mismatch_rate = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_MAX_MISMATCH_RATE, PindelNodeModel.DEF_MAX_MISMATCH_RATE, PindelNodeModel.MIN_MAX_MISMATCH_RATE, PindelNodeModel.MAX_MAX_MISMATCH_RATE);
	
	// pindel2vcf params
	final SettingsModelBoolean use_ref_filename = new SettingsModelBoolean(PindelNodeModel.CFGKEY_USE_REF_FILENAME, PindelNodeModel.DEF_USE_REF_FILENAME);
	final SettingsModelString refname = new SettingsModelString(PindelNodeModel.CFGKEY_REFNAME, PindelNodeModel.DEF_REFNAME);
	final SettingsModelBoolean use_cur_date = new SettingsModelBoolean(PindelNodeModel.CFGKEY_USE_CUR_DATE, PindelNodeModel.DEF_USE_CUR_DATE);
	final SettingsModelDate refdate = new SettingsModelDate(PindelNodeModel.CFGKEY_REFDATE);
	final SettingsModelIntegerBounded min_reads = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_READS, PindelNodeModel.DEF_MIN_READS, PindelNodeModel.MIN_MIN_READS, PindelNodeModel.MAX_MIN_READS);
	final SettingsModelDoubleBounded hetero_frac = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_HETERO_FRAC, PindelNodeModel.DEF_HETERO_FRAC, PindelNodeModel.MIN_HETERO_FRAC, PindelNodeModel.MAX_HETERO_FRAC);
	final SettingsModelDoubleBounded homo_frac = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_HOMO_FRAC, PindelNodeModel.DEF_HOMO_FRAC, PindelNodeModel.MIN_HOMO_FRAC, PindelNodeModel.MAX_HOMO_FRAC);
	final SettingsModelBoolean gatk_comp = new SettingsModelBoolean(PindelNodeModel.CFGKEY_GATK_COMP, PindelNodeModel.DEF_GATK_COMP);

	final SettingsModelIntegerBounded min_supp_reads = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_SUPP_READS, PindelNodeModel.DEF_MIN_SUPP_READS, PindelNodeModel.MIN_MIN_SUPP_READS, PindelNodeModel.MAX_MIN_SUPP_READS);
	final SettingsModelBoolean both_strands = new SettingsModelBoolean(PindelNodeModel.CFGKEY_BOTH_STRANDS, PindelNodeModel.DEF_BOTH_STRANDS);
	final SettingsModelIntegerBounded min_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_SIZE, PindelNodeModel.DEF_MIN_SIZE, PindelNodeModel.MIN_MIN_SIZE, PindelNodeModel.MAX_MIN_SIZE);
	final SettingsModelBoolean limit_size = new SettingsModelBoolean(PindelNodeModel.CFGKEY_LIMIT_SIZE, PindelNodeModel.DEF_LIMIT_SIZE);
	final SettingsModelIntegerBounded max_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MAX_SIZE, PindelNodeModel.DEF_MAX_SIZE, PindelNodeModel.MIN_MAX_SIZE, PindelNodeModel.MAX_MAX_SIZE);
	
	
    protected PindelNodeDialog() {
        super();
        
        GeneralOptions();
        PindelParams();
        Pindel2VCFParams();
    }
    
    private void Pindel2VCFParams(){
    	
    	createNewTab("Pindel2vcf Parameters");
    	
    	createNewGroup("Reference sequence");
    	addDialogComponent(new DialogComponentBoolean(use_ref_filename, "Use file name as reference name"));
    	addDialogComponent(new DialogComponentString(refname, "Name of the reference sequence:", true, 10));
    	refname.setEnabled(false);
    	
    	use_ref_filename.addChangeListener( new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				
				refname.setEnabled(!use_ref_filename.getBooleanValue());
				
			}
		});
    	
    	addDialogComponent(new DialogComponentBoolean(use_cur_date, "Use current date"));
    	addDialogComponent(new DialogComponentDate(refdate, "Date of the version of the reference sequence"));
    	refdate.setEnabled(false);
    	
    	use_cur_date.addChangeListener( new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				
				refdate.setEnabled(!use_cur_date.getBooleanValue());
				
			}
		});
    	
    	createNewGroup("Genotype");
    	addDialogComponent(new DialogComponentNumber(min_reads, "Minimum reads to report genotype", 1, 5));
    	addDialogComponent(new DialogComponentNumber(hetero_frac, "Proportion of reads defined as heterozygous", 0.01, 5));
    	addDialogComponent(new DialogComponentNumber(homo_frac, "Proportion of reads defined as homozygous", 0.01, 5));
    	addDialogComponent(new DialogComponentBoolean(gatk_comp, "Output GATK-compatible genotypes (recommended)"));
    	
    	createNewGroup("Filter");
    	addDialogComponent(new DialogComponentBoolean(both_strands, "Only output variants that are supported by reads on both strands"));
    	addDialogComponent(new DialogComponentNumber(min_supp_reads, "Minimum number of supporting reads", 1, 5));
    	addDialogComponent(new DialogComponentNumber(min_size, "Minimum variant size", 1, 5));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(limit_size, "Limit variant size"));
    	addDialogComponent(new DialogComponentNumber(max_size, "Maximum variant size", 1, 5));
    	max_size.setEnabled(false);
    	setHorizontalPlacement(false);
    	
    	limit_size.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				
				max_size.setEnabled(limit_size.getBooleanValue());
			}
		});
    }
    
    private void GeneralOptions(){
    	
        createNewGroup("Path to Pindel executable");
        addDialogComponent(new DialogComponentFileChooser(pindel, "pindel", JFileChooser.OPEN_DIALOG, false));
        
        createNewGroup("Interval for variant calling");
        addDialogComponent(new DialogComponentBoolean(interval, "Restrict variant calling to a certain genomic region"));
        setHorizontalPlacement(true);
        addDialogComponent(new DialogComponentString(chrom, "Chromosome"));
        addDialogComponent(new DialogComponentNumberEdit(start, "Start", 8));
        addDialogComponent(new DialogComponentNumberEdit(end, "End", 8));
        chrom.setEnabled(false);
        start.setEnabled(false);
        end.setEnabled(false);
        setHorizontalPlacement(false);
        addDialogComponent(new DialogComponentLabel("(Chromosome name has to match reference and BAM file header)"));
        
        // text/number fields for chromosome, start and end only active if interval is used
        interval.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				chrom.setEnabled(interval.getBooleanValue());
				start.setEnabled(interval.getBooleanValue());
				end.setEnabled(interval.getBooleanValue());
			}
		});
        
        createNewGroup("Path to Pindel config file");
        addDialogComponent(new DialogComponentFileChooser(config_file, "pindelconfig", JFileChooser.OPEN_DIALOG, false));
        create_config.setEnabled(false);
        addDialogComponent(new DialogComponentBoolean(create_config, "Create config file (requires PicardTools: CollectInsertMetrics as previous node)"));
        
        // disables file chooser if config file is create by this node
        create_config.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				config_file.setEnabled(!create_config.getBooleanValue());
			}
		});
        
        createNewGroup("Output");
        addDialogComponent(new DialogComponentBoolean(vcf_out, "Convert Pindel output (deletions and small insertions) to VCF format"));
        DialogComponentFileChooser p2vcf_fc= new DialogComponentFileChooser(pindel2vcf, "p2vcf", JFileChooser.OPEN_DIALOG);
        p2vcf_fc.setBorderTitle("Path to Pindel2vcf converter");
        addDialogComponent(p2vcf_fc);
        
        // disable file chooser for pindel2vcf converter if output is not converted
        vcf_out.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				pindel2vcf.setEnabled(vcf_out.getBooleanValue());
				setEnabled(vcf_out.getBooleanValue(), "Pindel2vcf parameters");
			}
		});
        
        createNewGroup("Runtime and memeroy usage");
        setHorizontalPlacement(true);
        addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1, 5));
        addDialogComponent(new DialogComponentNumber(bin_size, "Bin size", 1, 5));
        addDialogComponent(new DialogComponentLabel("(Mb of reference in memory)"));
        setHorizontalPlacement(false);
    	
    }
    
    private void PindelParams(){
    	
        createNewTab("Pindel Parameters");
        
        createNewGroup("Minimum number of matching bases");
        addDialogComponent(new DialogComponentLabel("Only consider reads as evidence if they map with more than this number of bases:"));
        addDialogComponent(new DialogComponentNumber(min_match_bases, "Mismatching bases per read", 1, 5));
        
        createNewGroup("Mismatch threshold");
        addDialogComponent(new DialogComponentLabel("Do not align a read if there is another mapping position below this threshold:"));
        addDialogComponent(new DialogComponentNumber(additional_mismatch, "Mismatching bases per alignment", 1, 5));
        
        createNewGroup("Number of perfect matches at breakpoints");
        addDialogComponent(new DialogComponentLabel("Number of perfectly matching bases around a breakpoint of a split read:"));
        addDialogComponent(new DialogComponentNumber(min_match_breakpoint, "Number of perfect matches", 1, 5));
        
        createNewGroup("Sequencing error rate");
        addDialogComponent(new DialogComponentLabel("Expected fraction of sequencing errors:"));
        addDialogComponent(new DialogComponentNumber(seq_err, "Error rate", 0.01, 5));
        createNewGroup("Maximum allowed mismatch rate");
        
        addDialogComponent(new DialogComponentLabel("Only consider aligned reads with mismatch rate below this fraction:"));
        addDialogComponent(new DialogComponentNumber(max_mismatch_rate, "Fraction of mismatching bases per read", 0.01, 5));
    	
    }
    
    
    /* 
     * Program:   pindel2vcf (conversion of Pindel output to VCF format)
Example:   pindel2vcf -p sample3chr20_D -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36
              -d 20101123-v sample3chr20_D.vcf

Note:      -is only guaranteed to work correctly on output files produced by pindel version 0.2.3 and above.
           -LI and BP files (long insertion and break point files) have a different type of header and
            are not supported yet.

-r/--reference  The name of the file containing the reference genome: required parameter
*-R/--reference_name  The name and version of the reference genome: required parameter
*-d/--reference_date  The date of the version of the reference genome used: required parameter
-p/--pindel_output  The name of the pindel output file containing the SVs
-P/--pindel_output_root  The root-name of the pindel output file; this will result in one big output file containing                                                     deletions, short and long insertions, tandem duplications and inversions
-v/--vcf  The name of the output vcf-file (default: name of pindel output file +".vcf"
-c/--chromosome  The name of the chromosome (default: SVs on all chromosomes are processed)
-w/--window_size  Memory saving option: the size of the genomic region in a chromosome of which structural variants are calculated separately, in millions of bases (default 300, for memory saving 100 or 50 recommended)
*-mc/--min_coverage  The minimum number of reads to provide a genotype (default 10)
*-he/--het_cutoff  The propertion of reads to call het (default 0.2)
*-ho/--hom_cutoff  The propertion of reads to call het (default 0.8)
*-is/--min_size  The minimum size of events to be reported (default 1)
*-as/--max_size  The maximum size of events to be reported (default infinite)
*-b/--both_strands_supported  Only report events that are detected on both strands (default false)
-m/--min_supporting_samples  The minimum number of samples an event needs to occur in with sufficient support to be reported (default 0)
*-e/--min_supporting_reads  The minimum number of supporting reads required for an event to be reported (default 1)
-f/--max_supporting_reads  The maximum number of supporting reads allowed for an event to be reported, allows protection against miscalls in due to segmental duplications or poorly mapped regions (default infinite)
-sr/--region_start  The start of the region of which events are to be reported (default 0)
-er/--region_end  The end of the region of which events are to be reported (default infinite)
-ir/--max_internal_repeats  Filters out all indels where the inserted/deleted sequence is a homopolymer/microsatellite of more than X repetitions (default infinite). For example: T->TCACACA has CACACA as insertion, which is a microsattelite of 3 repeats; this would be filtered out by setting -ir to 2
-co/--compact_output_limit  Puts all structural variations of which either the ref allele or the alt allele exceeds the specified size (say 10 in '-co 10') in the format 'chrom pos first_base <SVType>'
-il/--max_internal_repeatlength  Filters out all indels where the inserted/deleted sequence is a homopolymers/microsatellite with an unit size of more than Y, combine with the option -ir. Default value of -il is infinite. For example: T->TCAGCAG has CAGCAG as insertion, which has the fundamental repetitive unit CAG of length 3. This would be filtered out if -il has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -il is 2
-pr/--max_postindel_repeats  Filters out all indels where the inserted/deleted sequence is followed by a repetition (of over X times) of the fundamental repeat unit of the inserted/deleted sequence. For example, T->TCACA would usually be a normal insertion, which is not filtered out, but if the real sequence change is TCACACA->TCACACACACA, it will be filtered out by -pr of 1 or above, as the fundamental repeat unit of the inserted sequence (CA) is repeated more than one time in the postindel sequence [indel sequence CACA, postindel sequence CACACA]. Note: when CAC is inserted next to ACACAC, the repeat sequence is recognized as CA, even though the 'postrepeat' sequence is ACACAC
-pl/--max_postindel_repeatlength  Filters out all indels where the inserted/deleted sequence is followed by a repetition of  the fundamental repeat unit of the inserted/deleted sequence; the maximum size of that 'fundamental unit' given by the value of -pl (default infinite) For example: TCAG->TCAGCAG has insertion CAG and post-insertion sequence CAG. This insertion would be filtered out if -pl has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -pl is 2
-sb/--only_balanced_samples  Only count a sample as supporting an event if it is supported by reads on both strands, minimum reads per strand given by the -ss parameter. (default false)
-ss/--minimum_strand_support  Only count a sample as supporting an event if at least one of its strands is supported by X reads (default 1)
*-G/--gatk_compatible  calls genotypes which could either be homozygous or heterozygous not as ./1 but as 0/1, to ensure compatibility with GATK
     * 
     */
}

