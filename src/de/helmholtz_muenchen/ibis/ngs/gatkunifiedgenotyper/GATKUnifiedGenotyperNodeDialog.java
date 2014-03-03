package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GATKUnifiedGenotyper" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class GATKUnifiedGenotyperNodeDialog extends DefaultNodeSettingsPane {
	
	
	/* main options
	 * 
	 * path to gatk
	 * dbsnp file
	 * interval file
	 * SNP vs. INDEL -> genotype-likelihoods model
	 * #threads
	 */
	
	
	final SettingsModelString gatk = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_GATK, "");
	final SettingsModelBoolean use_dbsnp=new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USE_DBSNP, GATKUnifiedGenotyperNodeModel.DEF_USE_DBSNP);
	final SettingsModelString dbsnp_file=new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_DBSNP_FILE, GATKUnifiedGenotyperNodeModel.DEF_DBSNP_FILE);
	final SettingsModelBoolean use_interval = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USE_INTERVAL, GATKUnifiedGenotyperNodeModel.DEF_USE_INTERVAL);
	final SettingsModelString interval_file = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_INTERVAL_FILE, GATKUnifiedGenotyperNodeModel.DEF_INTERVAL_FILE);
	final SettingsModelString variant_type= new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_VARIANT_TYPE, GATKUnifiedGenotyperNodeModel.DEF_VARIANT_TYPE);
	final SettingsModelIntegerBounded num_threads = new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_NUM_THREADS, GATKUnifiedGenotyperNodeModel.DEF_NUM_THREADS, GATKUnifiedGenotyperNodeModel.MIN_NUM_THREADS, GATKUnifiedGenotyperNodeModel.MAX_NUM_THREADS);
	
	
	/* other options
	 * 
	 * !!!double stand_call_conf (30)
	 * !!!double stand_emit_conf (30)
	 *  
	 *  snps only
	 *  int minimum base quality required for calling (17)
	 *  
	 *  double pcr error rate (1.0e-4)
	 *  double contamination (0.0) -> tries to remove reads
	 *  double heterozygosity (0.001)
	 *  double max_deletion_fraction : more than this fraction reads with deletion -> base is not called (0.05)
	 *  
	 *  indels only
	 *  double indel_heterozygosity (1.25e-4)
	 *  int minindelcount (5)  
	 *  int minindelfrac (0.25)	 * 
	 *  byte indel gap open penalty (45)
	 *  byte indel gap extension penalty (10)
	 */
	
	final SettingsModelBoolean mbq = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_MBQ, GATKUnifiedGenotyperNodeModel.DEF_MBQ);
	final SettingsModelString baq = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_BAQ, GATKUnifiedGenotyperNodeModel.DEF_BAQ);
	final SettingsModelDoubleBounded call_min_confidence= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.DEF_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MIN_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MAX_CALL_MIN_CONFIDENCE);
	final SettingsModelDoubleBounded emit_min_confidence = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.DEF_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MIN_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MAX_EMIT_MIN_CONFIDENCE);
	final SettingsModelDoubleBounded pcr_error = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_PCR_ERR, GATKUnifiedGenotyperNodeModel.DEF_PCR_ERR, GATKUnifiedGenotyperNodeModel.MIN_PCR_ERR, GATKUnifiedGenotyperNodeModel.MAX_PCR_ERR);
	final SettingsModelDoubleBounded contamination = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_CONTAMINATION, GATKUnifiedGenotyperNodeModel.DEF_CONTAMINATION, GATKUnifiedGenotyperNodeModel.MIN_CONTAMINATION, GATKUnifiedGenotyperNodeModel.MAX_CONTAMINATION);
	final SettingsModelDoubleBounded heterozygosity= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_HET, GATKUnifiedGenotyperNodeModel.DEF_HET, GATKUnifiedGenotyperNodeModel.MIN_HET, GATKUnifiedGenotyperNodeModel.MAX_HET);
	final SettingsModelDoubleBounded max_deletion_fraction= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.DEF_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.MIN_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.MAX_MAX_DELETION_FRAC);
	final SettingsModelIntegerBounded min_base_qual= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.DEF_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.MIN_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.MAX_MIN_BASE_QUAL);
	final SettingsModelDoubleBounded indel_heterozygosity= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_INDEL_HET, GATKUnifiedGenotyperNodeModel.DEF_INDEL_HET, GATKUnifiedGenotyperNodeModel.MIN_INDEL_HET, GATKUnifiedGenotyperNodeModel.MAX_INDEL_HET);
	final SettingsModelIntegerBounded min_indel_count = new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.DEF_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.MIN_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.MAX_MIN_INDEL_CNT);
	final SettingsModelDoubleBounded min_indel_frac = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.DEF_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.MIN_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.MAX_MIN_INDEL_FRAC);
	final SettingsModelIntegerBounded gap_open_pen= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.DEF_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.MIN_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.MAX_GAP_OPEN_PEN);
	final SettingsModelIntegerBounded gap_cont_pen= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.DEF_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.MIN_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.MAX_GAP_CONT_PEN);
	
	
	private boolean dbsnp=false;
	
    /**
     * New pane for configuring GATKUnifiedGenotyper node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected GATKUnifiedGenotyperNodeDialog() {
        super();
        
        generalOptions();
        UGOptions();
     
    }
    
    private void UGOptions() {
    	
        createNewTab("UnifiedGenotyper");
        createNewGroup("Phred-scaled confidence threshold");
        addDialogComponent(new DialogComponentNumber(call_min_confidence, "Confidence threshold for calling", 1, 5));
        addDialogComponent(new DialogComponentNumber(emit_min_confidence, "Confidence threshold for emitting", 1, 5));
        
        createNewGroup("PCR error");
        addDialogComponent(new DialogComponentNumber(pcr_error, "Error rate", 0.0001, 5));
        
        createNewGroup("Sample contamination");
        addDialogComponent(new DialogComponentNumber(contamination, "Contamination fraction of reads", 0.0001, 5));
        
        createNewGroup("Heterozygosity");
        addDialogComponent(new DialogComponentNumber(heterozygosity, "Heterozygosity value", 0.0001, 5));
        
        createNewGroup("Fraction of deletions");
        addDialogComponent(new DialogComponentNumber(max_deletion_fraction, "Maximum fraction of deletions for locus to be callable", 0.001, 5));
         
        createNewGroup("SNP calling");
        addDialogComponent(new DialogComponentNumber(min_base_qual, "Minimum base quality score for calling", 1, 5));
        
        createNewGroup("Indel calling");
        addDialogComponent(new DialogComponentNumber(indel_heterozygosity, "Heterozygosity value",0.00001 , 5));
        addDialogComponent(new DialogComponentNumber(min_indel_count, "Minimum count of indel reads", 1, 5));
        addDialogComponent(new DialogComponentNumber(min_indel_frac, "Minimum fraction of indel reads ", 0.01, 5));
        addDialogComponent(new DialogComponentNumber(gap_open_pen, "Indel gap open penalty", 1, 5));
        addDialogComponent(new DialogComponentNumber(gap_cont_pen, "Indel continuation penalty", 1, 5));
        
        createNewGroup("Per-base alignment qualities (BAQ)");
        addDialogComponent(new DialogComponentStringSelection(baq, "Choose BAQ mode ", GATKUnifiedGenotyperNodeModel.AVAIL_BAQ));
        
        createNewGroup("Malformed read filter");
        addDialogComponent(new DialogComponentBoolean(mbq, "Filter reads with mismatching number of bases and base qualities"));
		
	}

	private void generalOptions(){
        // gatk executable
        createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(gatk, "gatk3", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);
    	
    	// dbsnp set for annotation
    	createNewGroup("SNPs from DBSNP");
        addDialogComponent(new DialogComponentBoolean(use_dbsnp, "Annotate SNPs with ID from DBSNP"));
        DialogComponentFileChooser dbsnpf=new DialogComponentFileChooser(dbsnp_file, "dbsnpset2", JFileChooser.OPEN_DIALOG, false, ".vcf");
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
        createNewGroup("Interval for variant calling");
        addDialogComponent(new DialogComponentBoolean(use_interval, "Restrict variant discovery to certain genomic regions"));
        interval_file.setEnabled(false);
        addDialogComponent(new DialogComponentFileChooser(interval_file, "ifile3", JFileChooser.OPEN_DIALOG, false, ".bed", ".intervals"));
        
        // add change listener that enables file chooser for interval file
        use_interval.addChangeListener(new ChangeListener(){
        	
        	public void stateChanged(ChangeEvent e){
        		
        		interval_file.setEnabled(use_interval.getBooleanValue());
        	}
        	
        });
        
        // variant types
        createNewGroup("Varaints");
        addDialogComponent(new DialogComponentLabel("Choose the variant types to be called"));
        addDialogComponent(new DialogComponentButtonGroup(variant_type, "", true, new String[] {"SNPs","INDELs","SNPs and Indels (separate files)"}, GATKUnifiedGenotyperNodeModel.AVAIL_VARIANT_TYPE));
        
        //#threads
        createNewGroup("Number of threads");
        addDialogComponent(new DialogComponentNumber(num_threads, "Threads", 1));
    }
}

