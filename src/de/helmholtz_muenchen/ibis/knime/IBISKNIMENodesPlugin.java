
package de.helmholtz_muenchen.ibis.knime;

import java.nio.file.Files;
import java.nio.file.Paths;

import javax.swing.JOptionPane;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.core.node.NodeLogger;
import org.osgi.framework.BundleContext;

import de.helmholtz_muenchen.ibis.knime.preferences.KNIMEPreferencePage;

/**
 * This is the OSGI bundle activator.
 * 
 * @author hastreiter
 */
public class IBISKNIMENodesPlugin extends AbstractUIPlugin {
	
	public static final String REF_GENOME = "ref_genome";
	public static final String REF_GENOME_DEFAULT = "";
	
	public static final String USE_HTE = "hte";
	public static final boolean HTE_DEFAULT = false;
	
	public static final String THRESHOLD = "threshold";
	public static final int THRESHOLD_DEFAULT = 1;
	
	public static final String DB_FILE = "db_file";
	public static final String DB_FILE_DEFAULT = "";
	
	public static final String NOTIFY = "notify";
	public static final boolean NOTIFY_DEFAULT = false;
	
	public static final String EMAIL_HOST = "host";
	public static final String EMAIL_HOST_DEFAULT = "outmail.helmholtz-muenchen.de";
	
	public static final String EMAIL_SENDER = "sender";
	public static final String EMAIL_SENDER_DEFAULT = "ibis.knime@helmholtz-muenchen.de";

	public static final String EMAIL_RECEIVER = "receiver";
	public static final String EMAIL_RECEIVER_DEFAULT = "";

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
        for(String s: KNIMEPreferencePage.TOOLS.keySet()) {
			String path = this.getToolPathPreference(s);
			if(Files.notExists(Paths.get(path)) && !path.equals("")) {
				this.setToolPathPreference(s, "");
				warn_message += "Path to "+s+" has been changed!"+ System.getProperty("line.separator"); 
			}
		}
        
        if(warn_message.length()>0) {
        	warn_message += "Please edit the paths in the KNIME4NGS preference page!"+System.getProperty("line.separator");
        }
        
        //check path to reference genome
        String ref_genome = this.getRefGenomePreference();
        if(Files.notExists(Paths.get(ref_genome)) && !ref_genome.equals("")) {
			this.setRefGenomePreference("");
			warn_message += "Path to the reference genome has been changed!"+ 
					System.getProperty("line.separator") + 
					"You can reset it in the KNIME4NGS preference page!" +
					System.getProperty("line.separator"); 
		}
        
        //check path to database
        String database = this.getDBFilePreference();
        if(Files.notExists(Paths.get(database)) && !database.equals("")) {
			this.setDBFilePreference("");
			this.setHTEPreference(false);
			warn_message += "Path to the HTE database has been changed!"+ 
					System.getProperty("line.separator") + 
					"You can reset it or create a new one in the KNIME4NGS preference page!" +
					System.getProperty("line.separator"); 
		}
        
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
		store.setDefault(REF_GENOME, REF_GENOME_DEFAULT);
		store.setDefault(USE_HTE, HTE_DEFAULT);
		store.setDefault(THRESHOLD, THRESHOLD_DEFAULT);
		store.setDefault(DB_FILE, DB_FILE_DEFAULT);
		store.setDefault(NOTIFY, NOTIFY_DEFAULT);
		store.setDefault(EMAIL_HOST, EMAIL_HOST_DEFAULT);
		store.setDefault(EMAIL_SENDER, EMAIL_SENDER_DEFAULT);
		store.setDefault(EMAIL_RECEIVER, EMAIL_RECEIVER_DEFAULT);
	}
	
	public String getRefGenomePreference() {
		return getPreferenceStore().getString(REF_GENOME);
	}
	
	public void setRefGenomePreference(String path) {
		getPreferenceStore().setValue(REF_GENOME, path);
	}
	
	public String getToolPathPreference(String tool) {
		return getPreferenceStore().getString(tool);
	}
	
	public void setToolPathPreference(String tool, String path) {
		getPreferenceStore().setValue(tool, path);
	}
	
	public boolean getHTEPreference() {
		return getPreferenceStore().getBoolean(USE_HTE);
	}
	
	public void setHTEPreference(boolean use) {
		getPreferenceStore().setValue(USE_HTE, use);
	}
	
	public String getThresholdPreference() {
		return getPreferenceStore().getString(THRESHOLD);
	}
	
	public void setThresholdPreference(String t) {
		getPreferenceStore().setValue(THRESHOLD, t);
	}
	
	public String getDBFilePreference() {
		return getPreferenceStore().getString(DB_FILE);
	}
	
	public void setDBFilePreference(String s) {
		getPreferenceStore().setValue(DB_FILE, s);
	}

	public void setNotifyPreference(boolean b){
		getPreferenceStore().setValue(NOTIFY, b);
	}
	
	public boolean getNotifyPreference() {
		return getPreferenceStore().getBoolean(NOTIFY);
	}
	
	public String getEmailPreference(String id) {
		if(id.equals(EMAIL_HOST)) {
			return getPreferenceStore().getString(EMAIL_HOST);
		}
		if(id.equals(EMAIL_SENDER)) {
			return getPreferenceStore().getString(EMAIL_SENDER);
		}
		if(id.equals(EMAIL_RECEIVER)) {
			return getPreferenceStore().getString(EMAIL_RECEIVER);
		}
		return "";
	}
	
	public void setEmailPreference(String id, String s) {
		getPreferenceStore().setValue(id, s);
	}
}
