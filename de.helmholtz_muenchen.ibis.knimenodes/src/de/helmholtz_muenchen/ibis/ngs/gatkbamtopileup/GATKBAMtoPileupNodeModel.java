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


package de.helmholtz_muenchen.ibis.ngs.gatkbamtopileup;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of GATKBAMtoPileup.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKBAMtoPileupNodeModel extends NodeModel {
    
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
//	public static final String CFGKEY_INFILE = "INFILE";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
//	public static final String CFGKEY_BED_FILE = "BED";
	public static final String CFGKEY_GATK_NT = "NT";
	public static final String CFGKEY_GATK_NCT = "NCT";
	
	
    private final SettingsModelString m_GATK 				= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_GATK_PATH, "---");
//    private final SettingsModelString m_INFILE 				= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_INFILE, "---");
    private final SettingsModelString m_REF_GENOME 			= new SettingsModelString(GATKBAMtoPileupNodeModel.CFGKEY_REF_GENOME, "---");
//    private final SettingsModelOptionalString m_BED_FILE 	= new SettingsModelOptionalString(GATKBAMtoPileupNodeModel.CFGKEY_BED_FILE, "---",false);
    private final SettingsModelIntegerBounded m_GATK_NT 	= new SettingsModelIntegerBounded(GATKBAMtoPileupNodeModel.CFGKEY_GATK_NT, 1, 1, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_GATK_NCT 	= new SettingsModelIntegerBounded(GATKBAMtoPileupNodeModel.CFGKEY_GATK_NCT, 1, 1, Integer.MAX_VALUE);

	//The Output Col Names
	public static final String OUT_COL1 = "PILEUP";
	
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKBAMtoPileupNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected GATKBAMtoPileupNodeModel() {
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String INFILE = inData[0].iterator().next().getCell(0).toString();
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java");
    	command.add("-jar "+m_GATK.getStringValue());
    	command.add("-T Pileup");
      	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-I "+INFILE);
    	
    	String OUTFILE = INFILE+".pileup";
    	command.add("-o "+OUTFILE);
    	
    	command.add("-nt "+m_GATK_NT.getIntValue());
    	command.add("-nct "+m_GATK_NCT.getIntValue());
//    	if(m_BED_FILE.isActive()){
//    		command.add("-L "+m_BED_FILE.getStringValue());
//    	}
    	
    	
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE)};
    	
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

        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
		   	 m_GATK.saveSettingsTo(settings);
//		   	 m_INFILE.saveSettingsTo(settings);
		   	 m_REF_GENOME.saveSettingsTo(settings);
//		   	 m_BED_FILE.saveSettingsTo(settings);
		   	 m_GATK_NT.saveSettingsTo(settings);
		   	 m_GATK_NCT.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_GATK.loadSettingsFrom(settings);
//		   	 m_INFILE.loadSettingsFrom(settings);
		   	 m_REF_GENOME.loadSettingsFrom(settings);
//		   	 m_BED_FILE.loadSettingsFrom(settings);
		   	 m_GATK_NT.loadSettingsFrom(settings);
		   	 m_GATK_NCT.loadSettingsFrom(settings);
    }	

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
		   	 m_GATK.validateSettings(settings);
//		   	 m_INFILE.validateSettings(settings);
		   	 m_REF_GENOME.validateSettings(settings);
//		   	 m_BED_FILE.validateSettings(settings);
		   	 m_GATK_NT.validateSettings(settings);
		   	 m_GATK_NCT.validateSettings(settings);
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

