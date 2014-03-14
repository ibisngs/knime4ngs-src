package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKBaseRecalibration" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class GATKBaseRecalibrationNodeDialog extends DefaultNodeSettingsPane {
		
	// general options
	/* 
	 * String path to gatk
	 * boolean use 1000G phase 1 indels
	 * String path to 1000G phase 1 indels 
	 * boolean use mills and 1000G gold standard
	 * String path to mills and 1000G gold standard
	 * boolean use dbsnp
	 * String path to dbsnp
	 * boolean use interval file
	 * String path to interval file
	 * 
	 * boolean create before/after plots -> BaseRecalibrator (2) -> AnalyzeCovariates
	 */
	
	final SettingsModelString gatk = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_GATK, "");
	final SettingsModelBoolean use_phase1_1000G = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USE_PHASE1_1000G, GATKBaseRecalibrationNodeModel.DEF_USE_PHASE1_1000G);
	final SettingsModelString phase1_1000G_file = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PHASE1_1000G_FILE, GATKBaseRecalibrationNodeModel.DEF_PHASE1_1000G_FILE);
	final SettingsModelBoolean use_mills_1000G = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USE_MILLS_1000G, GATKBaseRecalibrationNodeModel.DEF_USE_MILLS_1000G);
	final SettingsModelString mills_1000G_file = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_MILLS_1000G_FILE, GATKBaseRecalibrationNodeModel.DEF_MILLS_1000G_FILE);
	final SettingsModelBoolean use_dbsnp=new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USE_DBSNP, GATKBaseRecalibrationNodeModel.DEF_USE_DBSNP);
	final SettingsModelString dbsnp_file=new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_DBSNP_FILE, GATKBaseRecalibrationNodeModel.DEF_DBSNP_FILE);
	final SettingsModelBoolean use_interval = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USE_INTERVAL, GATKBaseRecalibrationNodeModel.DEF_USE_INTERVAL);
	final SettingsModelString interval_file = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_INTERVAL_FILE, GATKBaseRecalibrationNodeModel.DEF_INTERVAL_FILE);
	final SettingsModelBoolean create_plots = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_CREATE_PLOTS, GATKBaseRecalibrationNodeModel.DEF_CREATE_PLOTS);
	final SettingsModelIntegerBounded cpu_threads= new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_CPU_THREADS, GATKBaseRecalibrationNodeModel.DEF_CPU_THREADS, GATKBaseRecalibrationNodeModel.MIN_CPU_THREADS, GATKBaseRecalibrationNodeModel.MAX_CPU_THREADS);
	
	/* 
	 * BaseRecalibrator
	 * 
	 * for every covariate boolean value (--covariate + --no_standard_covs)
	 * ContextCovariate  (StandardCovariate) 
	 * CycleCovariate  (StandardCovariate) 
	 * QualityScoreCovariate  (RequiredCovariate) 
	 * ReadGroupCovariate  (RequiredCovariate)
	 * RepeatLengthCovariate   
	 * RepeatUnitCovariate   
	 * RepeatUnitAndLengthCovariate
	 * 
	 * byte low quality tail, all read ends with quality below this threshold are not considered (2)
	 * double basgrp Open Gap penalty (40.0)
	 * 
	 * byte deletion default quality
	 * byte insertion default quality
	 * int indel context size (3) [1,13]
	 * byte mismatch default quality
	 * int mismatches context size (2) [1,13]
	 * 
	 * int maximum cycle value (500)
	 * 
	 * boolean simplify output, erase all tags except read group
	 * 
	 */
	
	final SettingsModelBoolean context_cov= new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_CONTEXT_COV, GATKBaseRecalibrationNodeModel.DEF_CONTEXT_COV);
	final SettingsModelBoolean cycle_cov= new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_CYCLE_COV, GATKBaseRecalibrationNodeModel.DEF_CYCLE_COV);
	final SettingsModelBoolean rep_len_cov= new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_REP_LEN_COV, GATKBaseRecalibrationNodeModel.DEF_REP_LEN_COV);
	final SettingsModelBoolean rep_unit_cov= new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_REP_UNIT_COV, GATKBaseRecalibrationNodeModel.DEF_REP_UNIT_COV);
	final SettingsModelIntegerBounded low_qual_tail= new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_LOW_QUAL_TAIL, GATKBaseRecalibrationNodeModel.DEF_LOW_QUAL_TAIL, GATKBaseRecalibrationNodeModel.MIN_LOW_QUAL_TAIL, GATKBaseRecalibrationNodeModel.MAX_LOW_QUAL_TAIL);
	final SettingsModelDoubleBounded gap_open = new SettingsModelDoubleBounded(GATKBaseRecalibrationNodeModel.CFGKEY_GAP_OPEN, GATKBaseRecalibrationNodeModel.DEF_GAP_OPEN, GATKBaseRecalibrationNodeModel.MIN_GAP_OPEN, GATKBaseRecalibrationNodeModel.MAX_GAP_OPEN);
	final SettingsModelIntegerBounded deletion_def_qual= new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_DELETION_DEF_QUAL, GATKBaseRecalibrationNodeModel.DEF_DELETION_DEF_QUAL, GATKBaseRecalibrationNodeModel.MIN_DELETION_DEF_QUAL, GATKBaseRecalibrationNodeModel.MAX_DELETION_DEF_QUAL);
	final SettingsModelIntegerBounded insertion_def_qual = new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_INSERTION_DEF_QUAL, GATKBaseRecalibrationNodeModel.DEF_INSERTION_DEF_QUAL, GATKBaseRecalibrationNodeModel.MIN_INSERTION_DEF_QUAL, GATKBaseRecalibrationNodeModel.MAX_INSERTION_DEF_QUAL);
	final SettingsModelIntegerBounded indel_context_size = new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_INDEL_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.DEF_INDEL_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.MIN_INDEL_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.MAX_INDEL_CONTEXT_SIZE);
	final SettingsModelIntegerBounded mismatch_def_qual  = new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_MISMATCH_DEF_QUAL, GATKBaseRecalibrationNodeModel.DEF_MISMATCH_DEF_QUAL, GATKBaseRecalibrationNodeModel.MIN_MISMATCH_DEF_QUAL, GATKBaseRecalibrationNodeModel.MAX_MISMACTH_DEF_QUAL);
	final SettingsModelIntegerBounded mismatch_context_size= new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_MISMATCH_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.DEF_MISMATCH_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.MIN_MISMATCH_CONTEXT_SIZE, GATKBaseRecalibrationNodeModel.MAX_MISMATCH_CONTEXT_SIZE);
	final SettingsModelIntegerBounded max_cycles= new SettingsModelIntegerBounded(GATKBaseRecalibrationNodeModel.CFGKEY_MAX_CYCLES, GATKBaseRecalibrationNodeModel.DEF_MAX_CYCLES, GATKBaseRecalibrationNodeModel.MIN_MAX_CYCLES, GATKBaseRecalibrationNodeModel.MAX_MAX_CYCLES);
	final SettingsModelBoolean simplify_out=new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_SIMPLIFY_OUT, GATKBaseRecalibrationNodeModel.DEF_SIMPLIY_OUT);
	
	private boolean phase1;
	private boolean mills;
	private boolean dbsnp;
	
	//Proxy options
	private final SettingsModelBoolean useproxy = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXY, false);
	final SettingsModelString proxyhost = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYHOST, null);
	final SettingsModelString proxyport = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPORT, null);
	private final SettingsModelBoolean useproxyauth = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXYAUTH, false);
	final SettingsModelString proxyuser = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYUSER, null);
	final SettingsModelString proxypassword = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPASSWORD, null);
	
	
	
    protected GATKBaseRecalibrationNodeDialog() {
        super();
        
        generateGeneralOptions();
        generateBaseRecalOptions();
        generateProxyOptions();
        

    }
    
    private void generateGeneralOptions(){
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(gatk, "gatk2", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
        // sets of known indels from database for realignment
        createNewGroup("Sets of known polymorphisms (at least one set has to be chosen)");
        
        addDialogComponent(new DialogComponentBoolean(use_phase1_1000G, "Use 1000 genomes phase1 indel set (required for proper recalibration)"));
        DialogComponentFileChooser p1f =new DialogComponentFileChooser(phase1_1000G_file, "p11000g2", JFileChooser.OPEN_DIALOG, false, ".vcf");
        p1f.setBorderTitle("Choose File (disabled if file available from previous node)");
        addDialogComponent(p1f);
        
        // add change listener that enables file chooser for 1000g phase1 indels
        use_phase1_1000G.addChangeListener(new ChangeListener() {
			
			public void stateChanged(ChangeEvent e) {
				
				if(!phase1){
					if(!use_phase1_1000G.getBooleanValue() && !phase1_1000G_file.isEnabled()){
						phase1=true;
					}
					else{
						phase1_1000G_file.setEnabled(use_phase1_1000G.getBooleanValue());
					}
				}
			}	
		});
        
        addDialogComponent(new DialogComponentBoolean(use_mills_1000G, "Use Mills and 1000 genomes gold standard indel set (required for proper recalibration)"));
        DialogComponentFileChooser millsf=new DialogComponentFileChooser(mills_1000G_file, "m1000ggs2", JFileChooser.OPEN_DIALOG, false, ".vcf");
        millsf.setBorderTitle("Choose File (disabled if file available from previous node)");
        addDialogComponent(millsf);        
        
        use_mills_1000G.addChangeListener(new ChangeListener() {
			
			public void stateChanged(ChangeEvent e) {
				
				if(!mills){
					if(!use_mills_1000G.getBooleanValue() && !mills_1000G_file.isEnabled()){
						mills=true;
					}
					else{
						mills_1000G_file.setEnabled(use_mills_1000G.getBooleanValue());
					}
				}
				
			}
		});
        
        addDialogComponent(new DialogComponentBoolean(use_dbsnp, "Use dbSNP snp set (required for proper recalibration)"));
        DialogComponentFileChooser dbsnpf=new DialogComponentFileChooser(dbsnp_file, "dbsnpset", JFileChooser.OPEN_DIALOG, false, ".vcf");
        dbsnpf.setBorderTitle("Choose File (disabled if file available from previous node)");
        addDialogComponent(dbsnpf);
        
        use_dbsnp.addChangeListener(new ChangeListener() {
			
			public void stateChanged(ChangeEvent e) {
				
				if(!dbsnp){
					if(!use_dbsnp.getBooleanValue() && !dbsnp_file.isEnabled()){
						dbsnp=true;
					}
					else{
						dbsnp_file.setEnabled(use_dbsnp.getBooleanValue());
					}
				}
				
			}
		});
        
        // interval file for realignment
        createNewGroup("Interval for recalibration");
        addDialogComponent(new DialogComponentBoolean(use_interval, "Restrict recalibration to certain genomic regions"));
        interval_file.setEnabled(false);
        addDialogComponent(new DialogComponentFileChooser(interval_file, "ifile2", JFileChooser.OPEN_DIALOG, false, ".bed", ".intervals"));
        
        // add change listener that enables file chooser for interval file
        use_interval.addChangeListener(new ChangeListener(){
        	
        	public void stateChanged(ChangeEvent e){
        		
        		interval_file.setEnabled(use_interval.getBooleanValue());
        	}
        	
        });
        
        createNewGroup("Before after plots");
        addDialogComponent(new DialogComponentBoolean(create_plots, "Create before after plots for base qualities"));
        
        createNewGroup("Number of cpu threads");
        addDialogComponent(new DialogComponentNumber(cpu_threads, "Number of cpu threads", 1, 5));
        
    }
        
    private void generateBaseRecalOptions(){
        createNewTab("BaseReclibrator");
        
        createNewGroup("Calculation of covariates");
        addDialogComponent(new DialogComponentLabel("Quality score covariate"));
        addDialogComponent(new DialogComponentLabel("Read group covariate"));
        addDialogComponent(new DialogComponentBoolean(context_cov, "Context covariate"));
        addDialogComponent(new DialogComponentBoolean(cycle_cov, "Cycle covariate"));
        addDialogComponent(new DialogComponentBoolean(rep_len_cov, "Repeat length covariate"));
        addDialogComponent(new DialogComponentBoolean(rep_unit_cov, "Repeat unit covariate"));
        
        createNewGroup("Gap open penalty");
        addDialogComponent(new DialogComponentNumber(gap_open, "Gap open penalty", 0.1, 5));

        createNewGroup("Quality threshold for read tails");
        addDialogComponent(new DialogComponentNumber(low_qual_tail, "Quality threshold", 1, 5));
        
        createNewGroup("Deletion, insertion and mismatch settings");
        addDialogComponent(new DialogComponentNumber(deletion_def_qual, "Default quality for deletions (set -1 to disable)", 1, 5));
        addDialogComponent(new DialogComponentNumber(insertion_def_qual, "Default quality for insertions (set -1 to disable)", 1, 5));
        addDialogComponent(new DialogComponentNumber(mismatch_def_qual, "Default quality for mismatches (set -1 to disable)", 1, 5));
        addDialogComponent(new DialogComponentNumber(indel_context_size, "k-mer context size for indels", 1,5));
        addDialogComponent(new DialogComponentNumber(mismatch_context_size,"k-mer context size for mismatches" , 1, 5));
        
        createNewGroup("Cycle Covariate");
        addDialogComponent(new DialogComponentNumber(max_cycles,"Cycle Threshold" , 1, 5));
        
        createNewGroup("Simplify output");
        addDialogComponent(new DialogComponentBoolean(simplify_out, "Remove all additional tags from BAM file (except read group)"));
        
        // disables max value for cycles if cycle covariate is not calculated
        cycle_cov.addChangeListener(new ChangeListener() {
			
			public void stateChanged(ChangeEvent e) {
				max_cycles.setEnabled(cycle_cov.getBooleanValue());
			}
		});
        
        // disables covariate choice
        create_plots.addChangeListener(new ChangeListener() {
			
			public void stateChanged(ChangeEvent e) {
				context_cov.setEnabled(!create_plots.getBooleanValue());
				cycle_cov.setEnabled(!create_plots.getBooleanValue());
				rep_len_cov.setEnabled(!create_plots.getBooleanValue());
				rep_unit_cov.setEnabled(!create_plots.getBooleanValue());
			}
		});
    }
    
    private void generateProxyOptions(){
	  	createNewTab("Proxy options");
	  	createNewGroup("General");
	  	addDialogComponent(new DialogComponentBoolean(useproxy, "Enable proxy"));
	  	addDialogComponent(new DialogComponentString(proxyhost, "Proxy host"));
	  	addDialogComponent(new DialogComponentString(proxyport, "Proxy port"));
	  	createNewGroup("Authentication");
	  	addDialogComponent(new DialogComponentBoolean(useproxyauth, "Enable authentication"));
	  	addDialogComponent(new DialogComponentString(proxyuser, "Proxy username"));
	  	addDialogComponent(new DialogComponentString(proxypassword, "Proxy password"));
	
	  	
	  	useproxy.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
						if(useproxy.getBooleanValue()){
							proxyhost.setEnabled(true);
							proxyport.setEnabled(true);
							useproxyauth.setEnabled(true);
						}else{
							proxyhost.setEnabled(false);
							proxyport.setEnabled(false);
							proxyuser.setEnabled(false);
							proxypassword.setEnabled(false);
							useproxyauth.setEnabled(false);
						}
				}
			});
	  	
	  	useproxyauth.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
						if(useproxy.getBooleanValue() && useproxyauth.getBooleanValue()){
							proxyuser.setEnabled(true);
							proxypassword.setEnabled(true);
						}else{
							proxypassword.setEnabled(false);
							proxyuser.setEnabled(false);
						}
				}
			});
  }
    
}

