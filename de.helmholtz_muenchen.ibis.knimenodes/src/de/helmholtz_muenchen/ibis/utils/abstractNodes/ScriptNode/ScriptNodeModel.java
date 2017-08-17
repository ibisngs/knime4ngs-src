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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import java.io.File;
import java.io.IOException;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;

public abstract class ScriptNodeModel extends HTExecutorNodeModel {
	protected final String SCRIPT;
	public static final String SCRIPTS_SUBFOLDER = "scripts";
	
	
	protected ScriptNodeModel(int nrInDataPorts, int nrOutDataPorts, String script, int numIn) {
		super(nrInDataPorts, nrOutDataPorts, numIn);		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected ScriptNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes, String script, int numIn) {
		super(inPortTypes, outPortTypes, numIn);		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected String getScriptPath(){
		return(IO.getScriptPath() + SCRIPTS_SUBFOLDER + File.separatorChar);
	}
	
	protected void executeScript(final ExecutionContext exec, String[] environment, File lockFile) throws Exception {
		executeCommand(this.getCommand(), this.getOutfile(), exec, environment, lockFile, null, null, null, null, null);
	}
	
	protected abstract String[] getCommand();
	protected abstract String getOutfile();

	public String getSCRIPT() {
		return SCRIPT;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
	}
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec);
    }
}
