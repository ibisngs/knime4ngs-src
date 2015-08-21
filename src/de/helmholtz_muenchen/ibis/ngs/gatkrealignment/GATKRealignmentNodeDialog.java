package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;



/**
 * <code>NodeDialog</code> for the "GATKRealignment" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class GATKRealignmentNodeDialog extends DefaultNodeSettingsPane {
	
	//options tab
	/*
	 * String path to gatk
	 * boolean use 1000G phase1 indels set
	 * String path to known indels: 1000G phase 1 indels
	 * boolean use Mills and 1000G gold standard indels set
	 * String path to known indels: Mills and 1000G gold standards indels
	 * boolean use interval file
	 * String path to interval file in bed format, perform only realignment in that interval
	 *  int number of threads 
	 */
	
	final SettingsModelString gatk = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_GATK, "");
	final SettingsModelBoolean use_phase1_1000G = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USE_PHASE1_1000G, GATKRealignmentNodeModel.DEF_USE_PHASE1_1000G);
	final SettingsModelString phase1_1000G_file = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_PHASE1_1000G_FILE, GATKRealignmentNodeModel.DEF_PHASE1_1000G_FILE);
	final SettingsModelBoolean use_mills_1000G = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USE_MILLS_1000G, GATKRealignmentNodeModel.DEF_USE_MILLS_1000G);
	final SettingsModelString mills_1000G_file = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_MILLS_1000G_FILE, GATKRealignmentNodeModel.DEF_MILLS_1000G_FILE);
	final SettingsModelBoolean use_interval = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USE_INTERVAL, GATKRealignmentNodeModel.DEF_USE_INTERVAL);
	final SettingsModelString interval_file = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_INTERVAL_FILE, GATKRealignmentNodeModel.DEF_INTERVAL_FILE);
	final SettingsModelIntegerBounded num_threads = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_NUM_THREADS, GATKRealignmentNodeModel.DEF_NUM_THREADS, GATKRealignmentNodeModel.MIN_NUM_THREADS, GATKRealignmentNodeModel.MAX_NUM_THREADS);
	final SettingsModelIntegerBounded memory_usage = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_JAVAMEMORY, GATKRealignmentNodeModel.DEF_NUM_JAVAMEMORY, GATKRealignmentNodeModel.MIN_NUM_JAVAMEMORY, GATKRealignmentNodeModel.MAX_NUM_JAVAMEMORY);
	
	//tab TargetCreator
	/*
	 * int max interval size for longer intervals no realignment will be performed because it is too time consuming (500)
	 * int min reads at locus to do entropy calculation (4)
	 * double mismatch fraction of base misamtches to have high entropy, only needed if alignment was performed by an ungapped aligner (0.0, 0<...<1)
	 * int window size, two SNPs or regions of high entropy within the window are clusterd together (10) (>1)
	 */
	
	final SettingsModelIntegerBounded max_interval = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_INTERVAL, GATKRealignmentNodeModel.DEF_MAX_INTERVAL, GATKRealignmentNodeModel.MIN_MAX_INTERVAL, GATKRealignmentNodeModel.MAX_MAX_INTERVAL);
	final SettingsModelIntegerBounded min_reads = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MIN_READS, GATKRealignmentNodeModel.DEF_MIN_READS, GATKRealignmentNodeModel.MIN_MIN_READS, GATKRealignmentNodeModel.MAX_MIN_READS);
	final SettingsModelIntegerBounded window = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_WINDOW, GATKRealignmentNodeModel.DEF_WINDOW, GATKRealignmentNodeModel.MIN_WINDOW, GATKRealignmentNodeModel.MAX_WINDOW);
	final SettingsModelDoubleBounded mismatch = new SettingsModelDoubleBounded(GATKRealignmentNodeModel.CFGKEY_MISMATCH, GATKRealignmentNodeModel.DEF_MISMATCH, GATKRealignmentNodeModel.MIN_MISMATCH, GATKRealignmentNodeModel.MAX_MISMATCH);

	
	//tab IndelRealigner
	/*
	 * String consensus determination model, how to compute consensus -> USE_READS, KNOWNS_ONLY, USE_SW (USE_READS)
	 * double LOD threshold for cleaning, significance threshold of position for realignment (5.0) 
	 * 
	 * advanced parameters -> reduce time for realignment for high coverage data
	 * double entropy threshold, percentage of mismatches at a genomic position, defines high entropy (0.15) 
	 * int maximal number of consensuses to try (30)
	 * int max insert size for realignment of paired-end reads (3000)
	 * int maximum number of base shifts allow when realigning a read (200)
	 * int max number of reads used for determining for consensuses (120)
	 * int max reads for realignemnt, if there are more reads no realignment will be performed (20000)
	 * 
	 * boolean do not output original alignemnt tags
	 */
	
	final SettingsModelString consensus_model = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_CONSENSUS_MODEL, GATKRealignmentNodeModel.DEF_CONSENSUS_MODEL);
	final SettingsModelDoubleBounded lod_threshold = new SettingsModelDoubleBounded(GATKRealignmentNodeModel.CFGKEY_LOD_THRESHOLD, GATKRealignmentNodeModel.DEF_LOD_THRESHOLD, GATKRealignmentNodeModel.MIN_LOD_THRESHOLD, GATKRealignmentNodeModel.MAX_LOD_THRESHOLD);
	final SettingsModelDoubleBounded entropy = new SettingsModelDoubleBounded(GATKRealignmentNodeModel.CFGKEY_ENTROPY, GATKRealignmentNodeModel.DEF_ENTROPY, GATKRealignmentNodeModel.MIN_ENTROPY, GATKRealignmentNodeModel.MAX_ENTROPY);
	final SettingsModelIntegerBounded max_consensuses = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_CONSENSUSES, GATKRealignmentNodeModel.DEF_MAX_CONSENSUSES, GATKRealignmentNodeModel.MIN_MAX_CONSENSUSES, GATKRealignmentNodeModel.MAX_MAX_CONSENSUSES);
	final SettingsModelIntegerBounded max_isize = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_ISIZE, GATKRealignmentNodeModel.DEF_MAX_ISIZE, GATKRealignmentNodeModel.MIN_MAX_ISIZE, GATKRealignmentNodeModel.MAX_MAX_ISIZE);
	final SettingsModelIntegerBounded max_pos_move = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_POS_MOVE, GATKRealignmentNodeModel.DEF_MAX_POS_MOVE, GATKRealignmentNodeModel.MIN_MAX_POS_MOVE, GATKRealignmentNodeModel.MAX_MAX_POS_MOVE);
	final SettingsModelIntegerBounded max_reads_cons = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_READS_CONS, GATKRealignmentNodeModel.DEF_MAX_READS_CONS, GATKRealignmentNodeModel.MIN_MAX_READS_CONS, GATKRealignmentNodeModel.MAX_MAX_READS_CONS);
	final SettingsModelIntegerBounded max_reads_realign = new SettingsModelIntegerBounded(GATKRealignmentNodeModel.CFGKEY_MAX_READS_REALIGN, GATKRealignmentNodeModel.DEF_MAX_READS_REALIGN, GATKRealignmentNodeModel.MIN_MAX_READS_REALIGN, GATKRealignmentNodeModel.MAX_MAX_READS_REALIGN);
	final SettingsModelBoolean alignment_tag = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_ALIGNMENT_TAG, GATKRealignmentNodeModel.DEF_ALIGNMENT_TAG);
	
	//indicate that files are available form previous node
	private boolean phase1;
	private boolean mills;
	
	//Proxy options
	private final SettingsModelBoolean useproxy = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXY, false);
	final SettingsModelString proxyhost = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_PROXYHOST, null);
	final SettingsModelString proxyport = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_PROXYPORT, null);
	private final SettingsModelBoolean useproxyauth = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXYAUTH, false);
	final SettingsModelString proxyuser = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_PROXYUSER, null);
	final SettingsModelString proxypassword = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_PROXYPASSWORD, null);
	
	public final SettingsModelOptionalString m_opt_flags = new SettingsModelOptionalString(GATKRealignmentNodeModel.CFGKEY_OPT_FLAGS,"",false);


    protected GATKRealignmentNodeDialog() {
        super();
        
        createGeneralOptions();
        createTargetCreatorOptions();        
        createIndelRealignmentOptions();
        createProxyOptions();
                    
    }
    
    private void createGeneralOptions(){
    	
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(gatk, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	gatkf.setBorderTitle("Choose File (disabled if file available from previous node)");
    	addDialogComponent(gatkf);

        // sets of known indels from database for realignment
        createNewGroup("Sets of known indels");
        addDialogComponent(new DialogComponentBoolean(use_phase1_1000G, "Use 1000 genomes phase1 indel set (improves quality of realignment)"));
        DialogComponentFileChooser p1f =new DialogComponentFileChooser(phase1_1000G_file, "p11000g", JFileChooser.OPEN_DIALOG, false, ".vcf");
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
        
        addDialogComponent(new DialogComponentBoolean(use_mills_1000G, "Use Mills and 1000 genomes gold standard indel set (improves quality of realignment)"));        
        DialogComponentFileChooser millsf=new DialogComponentFileChooser(mills_1000G_file, "m1000ggs", JFileChooser.OPEN_DIALOG, false, ".vcf");
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
        
        // interval file for realignment
        createNewGroup("Interval for realignment");
        addDialogComponent(new DialogComponentBoolean(use_interval, "Restrict realignment to certain genomic regions"));
        interval_file.setEnabled(false);
        addDialogComponent(new DialogComponentFileChooser(interval_file, "ifile", JFileChooser.OPEN_DIALOG, false, ".bed", ".intervals"));
        
        // add change listener that enables file chooser for interval file
        use_interval.addChangeListener(new ChangeListener(){
        	
        	public void stateChanged(ChangeEvent e){
        		
        		interval_file.setEnabled(use_interval.getBooleanValue());
        	}
        	
        });
        
        //#threads
        createNewGroup("Number of threads");
        addDialogComponent(new DialogComponentNumber(num_threads, "Threads", 1));
        
        //#Memory
        createNewGroup("Java Memory");
        addDialogComponent(new DialogComponentNumber(memory_usage, "Java Memory (GB) per thread", 1));
    }
    
    private void createTargetCreatorOptions(){
    	
        createNewTab("TargetCreator Options");
        
        createNewGroup("Maximum length of realignment interval");
        addDialogComponent(new DialogComponentNumber(max_interval, "Maximum interval length", 1, 5));
        
        createNewGroup("Minimum number of reads for entropy calculation");
        addDialogComponent(new DialogComponentNumber(min_reads, "Minimum number of reads", 1, 5));
        
        createNewGroup("Fraction of mismatching base qualities (only for ungapped aligner)");
        addDialogComponent(new DialogComponentNumber(mismatch, "Fraction of base quality mismatches at a position to have high entropy", 0.01, 5));
        
        createNewGroup("Window size for clustering SNPs");
        addDialogComponent(new DialogComponentNumber(window, "Window size", 1, 5));
    }
    
    private void createIndelRealignmentOptions(){
    	
        createNewTab("Realignment Options");
        
        createNewGroup("Choose how to determine consensus sequence");
        addDialogComponent(new DialogComponentStringSelection(consensus_model, "Consensus determination model", GATKRealignmentNodeModel.VALUES_CONSENSUS_MODEL));
        
        createNewGroup("Significance threshold (LOD Threshold) for realignment");
        addDialogComponent(new DialogComponentNumber(lod_threshold, "LOD Threshold", 0.1, 5));
        
        createNewGroup("Entropy threshold");
        addDialogComponent(new DialogComponentNumber(entropy, "Entropy threshold", 0.01, 5));
        
        createNewGroup("Maximum number of consensus sequences to try");
        addDialogComponent(new DialogComponentNumber(max_consensuses, "Consensus threshold", 1, 5));
        
        createNewGroup("Maximum insert size for realignment of paired-end reads");
        addDialogComponent(new DialogComponentNumber(max_isize, "Insert size threshold", 10, 5));
        
        createNewGroup("Maximum distance of read to original alignment position");
        addDialogComponent(new DialogComponentNumber(max_pos_move, "Read shift threshold", 1, 5));
        
        createNewGroup("Maximum number of reads used for consensus calculation");
        addDialogComponent(new DialogComponentNumber(max_reads_cons, "Reads for consensus threshold", 1, 5));
        
        createNewGroup("Maximum number of reads for realignment");
        addDialogComponent(new DialogComponentNumber(max_reads_realign, "Reads for realignemnt threshold", 100, 5));
        
        createNewGroup("Output");
        addDialogComponent(new DialogComponentBoolean(alignment_tag, "Do not output original cigar string"));
    
        createNewGroup("Further options");
        addDialogComponent(new DialogComponentOptionalString(m_opt_flags,"Further flags"));

    }
    
    private void createProxyOptions(){
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

