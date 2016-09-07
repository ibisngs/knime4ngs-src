/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.abstractNodes;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.port.PortType;

/**
 * <code>SettingsStorageNodeModel</code> can be used to store Stuff 
 * for the SettingsModel
 * @author Michael Kluge
 */
public abstract class SettingsStorageNodeModel extends NodeModel {
	
	// stores the defined SettingsModels objects
	private final HashMap<String, SettingsModel> SETTINGS_MAP = new HashMap<String, SettingsModel>();
	private final ArrayList<SettingsModel> SETTINGS = new ArrayList<SettingsModel>();
	
	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 */
	protected SettingsStorageNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
	}
	
	protected SettingsStorageNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes) {
		super(inPortTypes, outPortTypes);
	}
	
	/**
	 * should be used to add the settingModels to the SETTINGS hash
	 * in order to save and load ModelSettings automatically
	 * MUST BE CALLED IN CONSTRUCTOR explicitly! 
	 */
//	public abstract void init();
	
	public void addSetting(SettingsModel model) {
		if(model != null) {
			SETTINGS.add(model);
			// get protected config name
			try {
				Method m = model.getClass().getDeclaredMethod("getConfigName");
				m.setAccessible(true);
				String configName = m.invoke(model, new Object[]{}).toString();
				
				// add it to the settings map
				SETTINGS_MAP.put(configName, model);
			}
			catch(Exception e) {
//				e.printStackTrace();
			}
		}
		else {
			throw new IllegalArgumentException("SettingsModel can not be null if it should be added!");
		}
	}
	
	/**
	 * gets a saved settings or null, if not there
	 * @param name
	 */
	public SettingsModel getSetting(String name) {
		return this.SETTINGS_MAP.get(name);
	}
	
	
	/**************************************** KNIME METHODS ****************************************/
	/***********************************************************************************************/
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// get all models and validate them
    	for(SettingsModel set : SETTINGS)
    		set.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// get all models and load them
    	for(SettingsModel set : SETTINGS)
    		set.loadSettingsFrom(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	// get all models and save them
    	for(SettingsModel set : SETTINGS)
    		set.saveSettingsTo(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }
}
