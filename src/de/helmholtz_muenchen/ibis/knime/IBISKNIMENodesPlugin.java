
package de.helmholtz_muenchen.ibis.knime;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.preference.PreferenceConverter;
import org.eclipse.swt.graphics.RGB;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.core.node.NodeLogger;
import org.osgi.framework.BundleContext;

/**
 * This is the OSGI bundle activator.
 * 
 * @author hastreiter
 */
public class IBISKNIMENodesPlugin extends AbstractUIPlugin {
    
	public static final String FORMAT_PREFERENCE = "formatOnSave";
	public static final String INDENT_PREFERENCE = "indentSpaces";
	public static final String COMMENT_COLOR_PREFERENCE = "commentColor";
	public static final String ERROR_COLOR_PREFERENCE = "errorColor";
	public static final String VALID_COLOR_PREFERENCE = "validColor";
	
	public static final String TOOL_DIR_PREFERENCE = "tooldir";
	public static final String TOOL_DIR_DEFAULT = "~/";
	
	public static final String USE_HTE = "hte";
	public static final String THRESHOLD = "threshold";
	public static final String DB_FILE = "db_file";
	public static final String NOTIFY = "notify";
	public static final String EMAIL = "email";
	
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

//        IPreferenceStore store = IBISKNIMENodesPlugin.getDefault()
//                .getPreferenceStore();

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
		store.setDefault(FORMAT_PREFERENCE, true);
		store.setDefault(INDENT_PREFERENCE, 5);
		PreferenceConverter.setDefault(store, COMMENT_COLOR_PREFERENCE,  new RGB(0, 200, 125));
		PreferenceConverter.setDefault(store, ERROR_COLOR_PREFERENCE, new RGB(255, 0, 0));
		PreferenceConverter.setDefault(store, VALID_COLOR_PREFERENCE,  new RGB(0, 0, 0));

		store.setDefault(TOOL_DIR_PREFERENCE, TOOL_DIR_DEFAULT);
		store.setDefault(USE_HTE, false);
		store.setDefault(THRESHOLD, "1");
		store.setDefault(DB_FILE, "~/");
		store.setDefault(NOTIFY, false);
		store.setDefault(EMAIL, "");
	}
	
	
	public String getToolDirPreference() {
		return getPreferenceStore().getString(TOOL_DIR_PREFERENCE);
	}
	
	public void setToolDirPreference(String path) {
		getPreferenceStore().setValue(TOOL_DIR_PREFERENCE, path);
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
	
	public String getEmailPreference() {
		return getPreferenceStore().getString(EMAIL);
	}
	
	public void setEmailPreference(String s) {
		getPreferenceStore().setValue(EMAIL, s);
	}
}
