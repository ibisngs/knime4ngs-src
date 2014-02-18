package de.helmholtz_muenchen.ibis.utils.abstractNodes;

import java.util.HashMap;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
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
		
		// create initial HashMaps
		MODEL_SETTINGS_STRING_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelString>());
	}

	// stores the defined ModelSettings
	private static final HashMap<String, String> MODEL_SETTINGS_STRING_VALUES = new HashMap<String, String>();
	private static final HashMap<Integer, HashMap<String, SettingsModelString>> MODEL_SETTINGS_STRING_OBJECTS = new HashMap<Integer, HashMap<String, SettingsModelString>>();

	
	/************************************* SettingsModelString  *************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelString
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelString(String key, String defaultValue) {
		MODEL_SETTINGS_STRING_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelString Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelString getSettingsModelString(String key) throws IllegalArgumentException {
		return getSettingsModelString(key, this);
	}
	
	/**
	 * Returns a SettingsModelString Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelString object
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelString getSettingsModelString(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(MODEL_SETTINGS_STRING_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!MODEL_SETTINGS_STRING_OBJECTS.containsKey(instanceKey))
				MODEL_SETTINGS_STRING_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelString>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelString> models = MODEL_SETTINGS_STRING_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelString(key, MODEL_SETTINGS_STRING_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelString.")); 
	}
	
	
	/**************************************** KNIME METHODS ****************************************/
	/***********************************************************************************************/
	
    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// get all string models and validate them
    	for(SettingsModelString stringSetting : MODEL_SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		stringSetting.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// get all string models and load them
    	for(SettingsModelString stringSetting : MODEL_SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		stringSetting.loadSettingsFrom(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	// get all string models and save them
    	for(SettingsModelString stringSetting : MODEL_SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		stringSetting.saveSettingsTo(settings);
    }
}
