package de.helmholtz_muenchen.ibis.utils.abstractNodes;

import java.util.HashMap;

import org.knime.core.node.NodeModel;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>SettingsStorageNodeModel</code> can be used to store Stuff 
 * for the SettingsModel
 * @author Michael Kluge
 * TODO: support more cases than SettingsModelString
 */
public abstract class SettingsStorageNodeModel extends NodeModel {

	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 */
	protected SettingsStorageNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
	}

	// stores the defined ModelSettingsString pairs
	private static final HashMap<String, String> MODEL_SETTINGS_STRING = new HashMap<String, String>();
	
	/**
	 * adds a new SettingsModelString
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelString(String key, String defaultValue) {
		MODEL_SETTINGS_STRING.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelString Object
	 * @param key key of the SettingsModelString object
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelString getSettingsModelString(String key) throws IllegalArgumentException {
		// check, if defined
		if(MODEL_SETTINGS_STRING.containsKey(key))
			return new SettingsModelString(key, MODEL_SETTINGS_STRING.get(key));
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelString.")); 
	}
}
