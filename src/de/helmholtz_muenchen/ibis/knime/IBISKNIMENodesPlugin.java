
package de.helmholtz_muenchen.ibis.knime;

import java.util.StringTokenizer;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.preference.PreferenceConverter;
import org.eclipse.swt.graphics.RGB;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.core.node.NodeLogger;
import org.osgi.framework.BundleContext;

//import com.genericworkflownodes.util.FileStashFactory;

//import de.helmholtz_muenchen.ibis.knime.preferences.PreferenceInitializer;

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
	public static final String FONT_PREFERENCE = "font";
	public static final String BROWSER_PREFERENCE = "webBrowser";
	
//	public static final String EXEMPT_TAGS_PREFERENCE = "exemptTags";
//	public static final String ERRORS_PREFERENCE = "errors";
	
	public static final String TOOL_DIR_PREFERENCE = "tooldir";
	public static final String TOOL_DIR_DEFAULT = "~/";

	
	

//	public static final String NO_ERRORS = "fieldEditors.noError";
//	public static final String ERROR_FOR_MISSING_CLOSING_TAG = "fieldEditors.errorForMissingClosingTag";
//
//	public static final String EXEMPT_TAGS_DEFAULT = "<P>;<BR>;<IMG>;";
	
	
	
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
        System.out.println("toggling Debug Mode");
    }

    /**
     * Sets the debug status of the plugin.
     * 
     * @param debugEnabled
     *            The new debug status.
     */
    public static void setDebug(final boolean debugEnabled) {
        IBISKNIMENodesPlugin.DEBUG = debugEnabled;
        System.out.println("setting Debug Mode :" + debugEnabled);
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

        IPreferenceStore store = IBISKNIMENodesPlugin.getDefault()
                .getPreferenceStore();
//        FileStashFactory.setTempParentDirectory(new File(store
//                .getString(PreferenceInitializer.PREF_FILE_STASH_LOCATION)));
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
		
		
//		store.setDefault(EXEMPT_TAGS_PREFERENCE, EXEMPT_TAGS_DEFAULT);
//		store.setDefault(ERRORS_PREFERENCE, NO_ERRORS);	

		store.setDefault(TOOL_DIR_PREFERENCE, TOOL_DIR_DEFAULT);

		
		
	}
	
//	public String[] getDefaultExemptTagsPreference() {
//		return convert(getPreferenceStore().getDefaultString(TOOL_DIR_PREFERENCE));
//	}

	public String getToolDirPreference() {
		return getPreferenceStore().getString(TOOL_DIR_PREFERENCE);
	}
	
//	private String[] convert(String preferenceValue) {
//		StringTokenizer tokenizer =
//			new StringTokenizer(preferenceValue, ",");
//		int tokenCount = tokenizer.countTokens();
//		String[] elements = new String[tokenCount];
//
//		for (int i = 0; i < tokenCount; i++) {
//			elements[i] = tokenizer.nextToken();
//		}
//
//		return elements;
//	}

	public void setToolDirPreference(String path) {
		getPreferenceStore().setValue(TOOL_DIR_PREFERENCE, path);
	}

}
