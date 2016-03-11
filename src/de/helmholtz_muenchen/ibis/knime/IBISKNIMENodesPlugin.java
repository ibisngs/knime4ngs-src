
package de.helmholtz_muenchen.ibis.knime;

import java.nio.file.Files;
import java.nio.file.Paths;

import javax.swing.JOptionPane;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.core.node.NodeLogger;
import org.osgi.framework.BundleContext;

/**
 * This is the OSGI bundle activator.
 * 
 * @author hastreiter
 */
public class IBISKNIMENodesPlugin extends AbstractUIPlugin {
	
//	fields of preference page
	public static final String USE_HTE = "hte";
	public static final String DB_FILE = "db_file";
	public static final String THRESHOLD = "threshold";
	public static final String NOTIFY = "notify";
	public static final String EMAIL_HOST = "host";
	public static final String EMAIL_SENDER = "sender";
	public static final String EMAIL_RECEIVER = "receiver";
	public static final String REF_GENOME = "reference sequence";
	public static final String RES_HAPMAP = "HapMap";
	public static final String RES_OMNI = "Omni";
	public static final String RES_1000G = "1000G";
	public static final String RES_DBSNP = "dbSNP";
	public static final String RES_MILLS = "Mills";
	
//	default values
	public static final boolean HTE_DEFAULT = false;
	public static final int THRESHOLD_DEFAULT = 1;
	public static final boolean NOTIFY_DEFAULT = false;
	public static final String EMAIL_HOST_DEFAULT = "outmail.helmholtz-muenchen.de";
	public static final String EMAIL_SENDER_DEFAULT = "ibis.knime@helmholtz-muenchen.de";
	public static final String EMAIL_RECEIVER_DEFAULT = "";
	
//	tool binaries
	public static final String BCFTOOLS = "bcftools";
	public static final String BOWTIE2 = "bowtie2";
	public static final String BWA = "bwa";
	public static final String FEATURE_COUNTS = "featureCounts";
	public static final String GATK = "GenomeAnalysisTK.jar";
	public static final String PINDEL = "pindel";
	public static final String PINDEL2VCF = "pindel2vcf";
	public static final String SAMTOOLS = "samtools";
	public static final String SEGEMEHL = "segemehl.x";
	public static final String STAR = "STAR";
	public static final String VCFTOOLS = "vcftools";
	public static final String VEP = "variant_effect_predictor.pl";
	public static final String VEP_FILTER = "filter_vep.pl";
	
//	tool dependencies, neither in FIELDS nor in PATHS
	public static final String BOWTIE2_BUILD = "bowtie2-build";
	
//	exluded binaries: integrate binaries by adding them to PATHS and TOOLS in the KNIMEPreferencePage
//	public static final String BFAST = "bfast";
	
	public static String [] FIELDS = {USE_HTE, THRESHOLD, DB_FILE, NOTIFY, EMAIL_HOST, EMAIL_SENDER, EMAIL_RECEIVER};
	public static String [] PATHS = {REF_GENOME, RES_HAPMAP, RES_OMNI, RES_1000G, RES_DBSNP, RES_MILLS, DB_FILE, BCFTOOLS, BOWTIE2, BWA, FEATURE_COUNTS, GATK, PINDEL, PINDEL2VCF, SAMTOOLS, SEGEMEHL, STAR, VCFTOOLS, VEP, VEP_FILTER};
	
	/**
     * The shared instance.
     */
    private static IBISKNIMENodesPlugin IKN_PLUGIN;

    /**
     * The central static logger.
     */
    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(IBISKNIMENodesPlugin.class);

    /**
     * Debuggin state of the plugin.
     */
    private static boolean DEBUG = false;

    /**
     * Logging method for debugging purpose.
     * 
     * @param message
     *            The message to log.
     */
    public static void log(final String message) {
        if (IBISKNIMENodesPlugin.DEBUG) {
            LOGGER.info(message);
        }
    }

    /**
     * Check if the plugin is in DEBUG mode.
     * 
     * @return True if debugging is enabled, false otherwise.
     */
    public static boolean isDebug() {
        return IBISKNIMENodesPlugin.DEBUG;
    }

    /**
     * Change debug setting.
     */
    public static void toggleDebug() {
        IBISKNIMENodesPlugin.DEBUG = !IBISKNIMENodesPlugin.DEBUG;
        log("toggling Debug Mode");
    }

    /**
     * Sets the debug status of the plugin.
     * 
     * @param debugEnabled
     *            The new debug status.
     */
    public static void setDebug(final boolean debugEnabled) {
        IBISKNIMENodesPlugin.DEBUG = debugEnabled;
        log("setting Debug Mode :" + debugEnabled);
    }

    /**
     * The constructor.
     */
    public IBISKNIMENodesPlugin() {
        super();
    }

    /**
     * This method is called upon plug-in activation.
     * 
     * @param context
     *            The OSGI bundle context
     * @throws Exception
     *             If this GKN_PLUGIN could not be started
     */
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        IKN_PLUGIN = this;

        log("starting IKN_PLUGIN: IBISKNIMENodesPlugin");

        String warn_message = "";
        
        //check tool paths
        for(String s: PATHS) {
			String path = getStringPreference(s);
			if(Files.notExists(Paths.get(path)) && !path.equals("")) {
				getDefault().getPreferenceStore().setToDefault(s);
				warn_message += "Path to "+s+" has been changed!"+ System.getProperty("line.separator"); 
				if(s.equals(DB_FILE)) {
					getDefault().getPreferenceStore().setToDefault(USE_HTE);
					getDefault().getPreferenceStore().setToDefault(NOTIFY);
					warn_message += "You can reset the path to hte.db or create a new database file in the KNIME4NGS preference page!" +
							System.getProperty("line.separator"); 
				}
			}
		}
        
        if(warn_message.length()>0) {
        	warn_message += "Please edit the paths in the KNIME4NGS preference page!"+System.getProperty("line.separator");
        }
        
        //check path to reference genome
//        String ref_genome = getStringPreference(REF_GENOME);
//        if(Files.notExists(Paths.get(ref_genome)) && !ref_genome.equals("")) {
//			getDefault().getPreferenceStore().setToDefault(REF_GENOME);
//			warn_message += "Path to the reference genome has been changed!"+ 
//					System.getProperty("line.separator") + 
//					"You can reset it in the KNIME4NGS preference page!" +
//					System.getProperty("line.separator"); 
//		}
        
        //check path to database
//        String database = getStringPreference(DB_FILE);
//        if(Files.notExists(Paths.get(database)) && !database.equals("")) {
//			getDefault().getPreferenceStore().setToDefault(DB_FILE);
//			
//			warn_message += "Path to the HTE database has been changed!"+ 
//					System.getProperty("line.separator") + 
//					"You can reset it or create a new one in the KNIME4NGS preference page!" +
//					System.getProperty("line.separator"); 
//		}
        
        if(warn_message.length()>0) {
        JOptionPane.showMessageDialog(null,
			    warn_message,
			    "Warning",
			    JOptionPane.WARNING_MESSAGE);
        }
    }

    /**
     * This method is called when the plug-in is stopped.
     * 
     * @param context
     *            The OSGI bundle context
     * @throws Exception
     *             If this GKN_PLUGIN could not be stopped
     */
    @Override
    public void stop(final BundleContext context) throws Exception {
    	IKN_PLUGIN = null;
        super.stop(context);
    }

    /**
     * Returns the shared instance.
     * 
     * @return Singleton instance of the Plugin
     */
    public static IBISKNIMENodesPlugin getDefault() {
        return IKN_PLUGIN;
    }

	/**
	 * Initializes a preference store with default preference values 
	 * for this plug-in.
	 * 
	 * @param store the preference store to fill
	 */
	protected void initializeDefaultPreferences(IPreferenceStore store) {
		store.setDefault(USE_HTE, HTE_DEFAULT);
		store.setDefault(THRESHOLD, THRESHOLD_DEFAULT);
		store.setDefault(NOTIFY, NOTIFY_DEFAULT);
		store.setDefault(EMAIL_HOST, EMAIL_HOST_DEFAULT);
		store.setDefault(EMAIL_SENDER, EMAIL_SENDER_DEFAULT);
		store.setDefault(EMAIL_RECEIVER, EMAIL_RECEIVER_DEFAULT);
		
		for(String t: PATHS) {
			store.setDefault(t, "");
		}
	}
	
//	public String getRefGenomePreference() {
//		return getPreferenceStore().getString(REF_GENOME);
//	}
	
//	public void setRefGenomePreference(String path) {
//		getPreferenceStore().setValue(REF_GENOME, path);
//	}
	
//	public String getToolPathPreference(String tool) {
//		return getPreferenceStore().getString(tool);
//	}
	
//	public void setToolPathPreference(String tool, String path) {
//		getPreferenceStore().setValue(tool, path);
//	}
	
//	public boolean getHTEPreference() {
//		return getPreferenceStore().getBoolean(USE_HTE);
//	}
//	
//	public void setHTEPreference(boolean use) {
//		getPreferenceStore().setValue(USE_HTE, use);
//	}
	
//	public String getThresholdPreference() {
//		return getPreferenceStore().getString(THRESHOLD);
//	}
	
//	public void setThresholdPreference(String t) {
//		getPreferenceStore().setValue(THRESHOLD, t);
//	}
	
//	public String getDBFilePreference() {
//		return getPreferenceStore().getString(DB_FILE);
//	}
	
//	public void setDBFilePreference(String s) {
//		getPreferenceStore().setValue(DB_FILE, s);
//	}

//	public void setNotifyPreference(boolean b){
//		getPreferenceStore().setValue(NOTIFY, b);
//	}
//	
//	public boolean getNotifyPreference() {
//		return getPreferenceStore().getBoolean(NOTIFY);
//	}
	
//	public String getEmailPreference(String id) {
//		if(id.equals(EMAIL_HOST)) {
//			return getPreferenceStore().getString(EMAIL_HOST);
//		}
//		if(id.equals(EMAIL_SENDER)) {
//			return getPreferenceStore().getString(EMAIL_SENDER);
//		}
//		if(id.equals(EMAIL_RECEIVER)) {
//			return getPreferenceStore().getString(EMAIL_RECEIVER);
//		}
//		return "";
//	}
	
//	public void setEmailPreference(String id, String s) {
//		getPreferenceStore().setValue(id, s);
//	}
	
	public static boolean getBooleanPreference(String id) {
		return IBISKNIMENodesPlugin.getDefault().getPreferenceStore().getBoolean(id);
	}
	
	public static void setBooleanPreference(String id, boolean value) {
		IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setValue(id, value);
	}
	
	public static String getStringPreference(String id) {
		return IBISKNIMENodesPlugin.getDefault().getPreferenceStore().getString(id);
	}
	
	public static void setStringPreference(String id, String value) {
		IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setValue(id, value);
	}
	
	public static void setAllFieldsToDefault() {
		for(String f: FIELDS) {
			IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setToDefault(f);
		}
		for(String s: PATHS) {
			IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setToDefault(s);
		}
	}
}
