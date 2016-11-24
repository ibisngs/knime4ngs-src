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

package de.helmholtz_muenchen.ibis.ngs.plotgatkgenotypeconcordance;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of PlotGATKGenotypeConcordance.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class PlotGATKGenotypeConcordanceNodeModel extends NodeModel {
    
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(PlotGATKGenotypeConcordanceNodeModel.class);
	//The Output Col Names
	public static final String OUT_COL1 = "SummaryPlots";
    protected PlotGATKGenotypeConcordanceNodeModel() {
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

        String INFILE = inData[0].iterator().next().getCell(0).toString();
        StringBuffer sysErr = new StringBuffer(50);
    	StringBuffer sysOut = new StringBuffer(50);
        
    	String EVAL = getAvailableFlowVariables().get("EVAL").getStringValue();
    	String COMP = getAvailableFlowVariables().get("COMP").getStringValue();
    	
    	
        /**
         * Split Output before plotting
         */
        ArrayList<String> command = new ArrayList<String>();
    	String path2splitScript = IO.getScriptPath() + "scripts/perl/splitGenotypeConcordanceOutput.pl";
    	command.add("perl");
    	command.add(path2splitScript);
    	command.add(INFILE);
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,sysOut,sysErr);
    	
    	/**
    	 * Plot Summary 
    	 */
    	
        command = new ArrayList<String>();
    	String path2PlotScript = IO.getScriptPath() + "scripts/R/plotting/plotGenotypeConcordanceOutput.R";
    	command.add("Rscript");
    	command.add(path2PlotScript);
    	command.add(INFILE);
    	command.add(EVAL.substring(EVAL.lastIndexOf("/")+1));
    	command.add(COMP.substring(COMP.lastIndexOf("/")+1));
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,sysOut,sysErr);
    	
    	String Outfile = INFILE+"_SummaryPlots.pdf";
    	
    	
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(Outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	if(!(getAvailableFlowVariables().containsKey("EVAL") && getAvailableFlowVariables().containsKey("COMP"))){
    		throw new InvalidSettingsException("Node has to be connected to GATKGenotypeConcordance.");
    	}
    	
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
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

}

