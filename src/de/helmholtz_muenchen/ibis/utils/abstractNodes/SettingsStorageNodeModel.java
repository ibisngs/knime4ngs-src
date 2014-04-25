package de.helmholtz_muenchen.ibis.utils.abstractNodes;

import java.util.HashMap;

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleRange;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelLong;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;
import org.knime.core.node.port.PortType;
import org.knime.core.util.Pair;

/**
 * <code>SettingsStorageNodeModel</code> can be used to store Stuff 
 * for the SettingsModel
 * @author Michael Kluge
 */
public abstract class SettingsStorageNodeModel extends NodeModel {
	
	// stores the defined SettingsModelBoolean objects
	private static final HashMap<String, Boolean> SETTINGS_BOOLEAN_VALUES = new HashMap<String, Boolean>();
	private static final HashMap<Integer, HashMap<String, SettingsModelBoolean>> SETTINGS_BOOLEAN_OBJECTS = new HashMap<Integer, HashMap<String, SettingsModelBoolean>>();
	
	// stores the defined SettingsModelString objects
	private static final HashMap<String, String> SETTINGS_STRING_VALUES = new HashMap<String, String>();
	private static final HashMap<String, String> SETTINGS_OPTIONAL_STRING_VALUES = new HashMap<String, String>();
	private static final HashMap<String, String[]> SETTINGS_STRING_ARRAY_VALUES = new HashMap<String, String[]>();
	private static final HashMap<Integer, HashMap<String, SettingsModelString>> SETTINGS_STRING_OBJECTS 				= new HashMap<Integer, HashMap<String, SettingsModelString>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelOptionalString>> SETTINGS_OPTIONAL_STRING_OBJECT	= new HashMap<Integer, HashMap<String, SettingsModelOptionalString>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelStringArray>> SETTINGS_STRING_ARRAY_OBJECTS 		= new HashMap<Integer, HashMap<String, SettingsModelStringArray>>();

	// stores the defined SettingsModelNumber objects
	private static final HashMap<String, Integer> SETTINGS_INTEGER_VALUES 		= new HashMap<String, Integer>();
	private static final HashMap<String, Double> SETTINGS_DOUBLE_VALUES 		= new HashMap<String, Double>();
	private static final HashMap<String, Long> SETTINGS_LONG_VALUES 			= new HashMap<String, Long>();
	private static final HashMap<String, Pair<Double, Double>> SETTINGS_DOUBLE_RANGE_VALUES 								= new HashMap<String, Pair<Double, Double>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelInteger>> SETTINGS_INTEGER_OBJECTS 					= new HashMap<Integer, HashMap<String, SettingsModelInteger>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelDouble>> SETTINGS_DOUBLE_OBJECTS 					= new HashMap<Integer, HashMap<String, SettingsModelDouble>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelLong>> SETTINGS_LONG_OBJECTS 						= new HashMap<Integer, HashMap<String, SettingsModelLong>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelDoubleRange>> SETTINGS_DOUBLE_RANGE_OBJECTS 			= new HashMap<Integer, HashMap<String, SettingsModelDoubleRange>>();
	private static final HashMap<Integer, HashMap<String, SettingsModelOptionalString>> SETTINGS_OPTIONAL_STRING_OBJECTS 	= new HashMap<Integer, HashMap<String, SettingsModelOptionalString>>();

	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 */
	protected SettingsStorageNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
		init();
	}
	
	protected SettingsStorageNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes) {
		super(inPortTypes, outPortTypes);
		init();
	}
	
	private void init(){
		// create initial HashMaps
		SETTINGS_BOOLEAN_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelBoolean>());
		SETTINGS_STRING_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelString>());
		SETTINGS_OPTIONAL_STRING_OBJECT.put(hashCode(), new HashMap<String, SettingsModelOptionalString>());
		SETTINGS_STRING_ARRAY_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelStringArray>());
		SETTINGS_INTEGER_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelInteger>());
		SETTINGS_DOUBLE_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelDouble>());
		SETTINGS_LONG_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelLong>());
		SETTINGS_DOUBLE_RANGE_OBJECTS.put(hashCode(), new HashMap<String, SettingsModelDoubleRange>());
		
        // create all the settings because otherwise they are not saved, if not used in the model in configure...
        for(String key : SETTINGS_BOOLEAN_VALUES.keySet())
			getSettingsModelBoolean(key);
        for(String key : SETTINGS_STRING_VALUES.keySet())
			getSettingsModelString(key);
        for(String key : SETTINGS_OPTIONAL_STRING_VALUES.keySet())
			getSettingsModelOptionalString(key);
        for(String key : SETTINGS_STRING_ARRAY_VALUES.keySet())
			getSettingsModelStringArray(key);
        for(String key : SETTINGS_INTEGER_VALUES.keySet())
			getSettingsModelInteger(key);
        for(String key : SETTINGS_DOUBLE_VALUES.keySet())
			getSettingsModelDouble(key);
        for(String key : SETTINGS_LONG_VALUES.keySet())
			getSettingsModelLong(key);
        for(String key : SETTINGS_DOUBLE_RANGE_VALUES.keySet())
			getSettingsModelDoubleRange(key);
	}
	/************************************* SettingsModelBoolean  *************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelBoolean
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelBoolean(String key, boolean defaultValue) {
		SETTINGS_BOOLEAN_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a getSettingsModelBoolean Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelBoolean getSettingsModelBoolean(String key) throws IllegalArgumentException {
		return getSettingsModelBoolean(key, this);
	}
	
	/**
	 * Returns a getSettingsModelBoolean Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelString object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelBoolean getSettingsModelBoolean(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_BOOLEAN_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();

			// create a new hashMap, if not there yet
			if(!SETTINGS_BOOLEAN_OBJECTS.containsKey(instanceKey))
				SETTINGS_BOOLEAN_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelBoolean>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelBoolean> models = SETTINGS_BOOLEAN_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelBoolean(key, SETTINGS_BOOLEAN_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelBoolean.")); 
	}
	


	/************************************* SettingsModelString  *************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelString
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelString(String key, String defaultValue) {
		SETTINGS_STRING_VALUES.put(key, defaultValue);
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
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelString getSettingsModelString(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_STRING_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_STRING_OBJECTS.containsKey(instanceKey))
				SETTINGS_STRING_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelString>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelString> models = SETTINGS_STRING_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelString(key, SETTINGS_STRING_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelString.")); 
	}
	
	
	/************************************* SettingsModelString  *************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelOptionalString
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelOptionalString(String key, String defaultValue) {
		SETTINGS_OPTIONAL_STRING_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelOptionalString Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelOptionalString getSettingsModelOptionalString(String key) throws IllegalArgumentException {
		return getSettingsModelOptionalString(key, this);
	}
	
	/**
	 * Returns a SettingsModelOptionalString Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelOptionalString object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelOptionalString getSettingsModelOptionalString(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_OPTIONAL_STRING_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_OPTIONAL_STRING_OBJECTS.containsKey(instanceKey))
				SETTINGS_OPTIONAL_STRING_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelOptionalString>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelOptionalString> models = SETTINGS_OPTIONAL_STRING_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelOptionalString(key, SETTINGS_OPTIONAL_STRING_VALUES.get(key), false));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelOptionalString.")); 
	}
	
	/********************************** SettingsModelStringArray  **********************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new addSettingsModelStringArray
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelStringArray(String key, String[] defaultValue) {
		SETTINGS_STRING_ARRAY_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelStringArray Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelStringArray getSettingsModelStringArray(String key) throws IllegalArgumentException {
		return getSettingsModelStringArray(key, this);
	}
	
	/**
	 * Returns a SettingsModelStringArray Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelStringArray object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelStringArray getSettingsModelStringArray(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_STRING_ARRAY_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_STRING_ARRAY_OBJECTS.containsKey(instanceKey))
				SETTINGS_STRING_ARRAY_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelStringArray>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelStringArray> models = SETTINGS_STRING_ARRAY_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelStringArray(key, SETTINGS_STRING_ARRAY_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelStringArray.")); 
	}
	
	/************************************* SettingsModelInteger  ************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelInteger
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelInteger(String key, int defaultValue) {
		SETTINGS_INTEGER_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelInteger Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelInteger getSettingsModelInteger(String key) throws IllegalArgumentException {
		return getSettingsModelInteger(key, this);
	}
	
	/**
	 * Returns a SettingsModelInteger Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelInteger object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelInteger getSettingsModelInteger(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_INTEGER_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_INTEGER_OBJECTS.containsKey(instanceKey))
				SETTINGS_INTEGER_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelInteger>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelInteger> models = SETTINGS_INTEGER_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelInteger(key, SETTINGS_INTEGER_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelInteger.")); 
	}
	
	/************************************* SettingsModelDouble  ************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelDouble
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelDouble(String key, double defaultValue) {
		SETTINGS_DOUBLE_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelDouble Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelDouble getSettingsModelDouble(String key) throws IllegalArgumentException {
		return getSettingsModelDouble(key, this);
	}
	
	/**
	 * Returns a SettingsModelDouble Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelDouble object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelDouble getSettingsModelDouble(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_DOUBLE_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_DOUBLE_OBJECTS.containsKey(instanceKey))
				SETTINGS_DOUBLE_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelDouble>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelDouble> models = SETTINGS_DOUBLE_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelDouble(key, SETTINGS_DOUBLE_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelDouble.")); 
	}
	
	/************************************* SettingsModelLong  ************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelLong
	 * @param key key for the setting
	 * @param defaultValue default value of the setting
	 */
	protected static void addSettingsModelLong(String key, long defaultValue) {
		SETTINGS_LONG_VALUES.put(key, defaultValue);
	}
	
	/**
	 * Returns a SettingsModelLong Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelLong getSettingsModelLong(String key) throws IllegalArgumentException {
		return getSettingsModelLong(key, this);
	}
	
	/**
	 * Returns a SettingsModelLong Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelLong object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelLong getSettingsModelLong(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_LONG_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_LONG_OBJECTS.containsKey(instanceKey))
				SETTINGS_LONG_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelLong>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelLong> models = SETTINGS_LONG_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) 
				models.put(key, new SettingsModelLong(key, SETTINGS_LONG_VALUES.get(key)));
			
			// return the object
			return models.get(key);
		}
		else
			throw(new IllegalArgumentException("Key '" + key + "' is not defined for SettingsModelLong.")); 
	}
	
	/************************************* SettingsModelString  *************************************/
	/***********************************************************************************************/
	
	/**
	 * adds a new SettingsModelDoubleRange
	 * @param key key for the setting
	 * @param minValue min value of the range
	 * @param minValue max value of the range
	 */
	protected static void addSettingsModelRange(String key, double minValue, double maxValue) {
		SETTINGS_DOUBLE_RANGE_VALUES.put(key, new Pair<Double, Double>(minValue, maxValue));
	}
	
	/**
	 * Returns a SettingsModelDoubleRange Object for the model node
	 * @param key
	 * @return
	 * @throws IllegalArgumentException
	 */
	public SettingsModelDoubleRange getSettingsModelDoubleRange(String key) throws IllegalArgumentException {
		return getSettingsModelDoubleRange(key, this);
	}
	
	/**
	 * Returns a SettingsModelDoubleRange Object 
	 * Should be used called if call comes not from the Model!
	 * @param key key of the SettingsModelDoubleRange object
	 * @param caller object, which called wants to get the settings
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static SettingsModelDoubleRange getSettingsModelDoubleRange(String key, final Object caller) throws IllegalArgumentException {
		// check, if defined
		if(SETTINGS_DOUBLE_RANGE_VALUES.containsKey(key)) {
			Integer instanceKey = caller.hashCode();
			
			// create a new hashMap, if not there yet
			if(!SETTINGS_DOUBLE_RANGE_OBJECTS.containsKey(instanceKey))
				SETTINGS_DOUBLE_RANGE_OBJECTS.put(instanceKey, new HashMap<String, SettingsModelDoubleRange>());
				
			// check, if settings model is there
			HashMap<String, SettingsModelDoubleRange> models = SETTINGS_DOUBLE_RANGE_OBJECTS.get(instanceKey);
			if(!models.containsKey(key)) {
				Pair<Double, Double> p = SETTINGS_DOUBLE_RANGE_VALUES.get(key);
				models.put(key, new SettingsModelDoubleRange(key, p.getFirst(), p.getSecond()));
			}
			
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
    	// get all models and validate them
    	for(SettingsModelBoolean set : SETTINGS_BOOLEAN_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelString set : SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelString set : SETTINGS_OPTIONAL_STRING_OBJECT.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelStringArray set : SETTINGS_STRING_ARRAY_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelInteger set : SETTINGS_INTEGER_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelDouble set : SETTINGS_DOUBLE_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelLong set : SETTINGS_LONG_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    	for(SettingsModelDoubleRange set : SETTINGS_DOUBLE_RANGE_OBJECTS.get(hashCode()).values())
    		set.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// get all models and load them
    	for(SettingsModelBoolean set : SETTINGS_BOOLEAN_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelString set : SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelString set : SETTINGS_OPTIONAL_STRING_OBJECT.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelStringArray set : SETTINGS_STRING_ARRAY_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelInteger set : SETTINGS_INTEGER_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelDouble set : SETTINGS_DOUBLE_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelLong set : SETTINGS_LONG_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    	for(SettingsModelDoubleRange set : SETTINGS_DOUBLE_RANGE_OBJECTS.get(hashCode()).values())
    		set.loadSettingsFrom(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	// get all models and save them
    	for(SettingsModelBoolean set : SETTINGS_BOOLEAN_OBJECTS.get(hashCode()).values()) {
    		set.saveSettingsTo(settings); 
    		System.out.println(set);
    	}
    	for(SettingsModelString set : SETTINGS_STRING_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelString set : SETTINGS_OPTIONAL_STRING_OBJECT.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelStringArray set : SETTINGS_STRING_ARRAY_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelInteger set : SETTINGS_INTEGER_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelDouble set : SETTINGS_DOUBLE_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelLong set : SETTINGS_LONG_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    	for(SettingsModelDoubleRange set : SETTINGS_DOUBLE_RANGE_OBJECTS.get(hashCode()).values())
    		set.saveSettingsTo(settings);
    }
}
