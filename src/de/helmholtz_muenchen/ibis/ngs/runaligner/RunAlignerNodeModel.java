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

package de.helmholtz_muenchen.ibis.ngs.runaligner;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataCell;
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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * This is the model implementation of RunAligner.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiter
 */
public class RunAlignerNodeModel extends NodeModel {
    
	
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(RunAlignerNodeModel.class);
	
	// name of the output variables
	public static final String OUT_COL1 = "Path2ReadFile1";
	public static final String OUT_COL2 = "Path2ReadFile2";
	
	public static final String CFGKEY_READSEQFILE = "readseqfile";
	public static final String CFGKEY_READSEQFILE2 = "readseqfile2";
	public static final String CFGKEY_READTYPE = "readType";
	
	private final SettingsModelString m_readseqfile = new SettingsModelString(RunAlignerNodeModel.CFGKEY_READSEQFILE,"");
	private final SettingsModelString m_readseqfile2 = new SettingsModelString(RunAlignerNodeModel.CFGKEY_READSEQFILE2,"");
	private final SettingsModelString m_readType = new SettingsModelString(RunAlignerNodeModel.CFGKEY_READTYPE, "single-end");
		
    /**
     * Constructor for the node model.
     */
    protected RunAlignerNodeModel() {
    	
        super(0, 1);
        
        m_readseqfile2.setEnabled(false);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
      	
    	/**Init and Get Values**/
    	String readFile1 = m_readseqfile.getStringValue();
    	String readFile2 = "";
    	String readType = m_readType.getStringValue();
    	if(m_readseqfile2.isEnabled() && readType.equals("paired-end")) {
    		readFile2 = m_readseqfile2.getStringValue();
    	}
    	
    	LOGGER.info("Reads File 1: " + readFile1);
    	LOGGER.info("Reads File 2: " + readFile2);
    	LOGGER.info("Read Type: " + readType);

		if(new File(readFile1).exists() && new File(readFile1).exists()){
			LOGGER.info("Sequence file(s) found.");
		}else{
			throw new IOException("One or both sequence files cannot be found!");
		}
    	
    	/**Push FlowVariables**/
    	String isBAM = "";
    	if(readFile1.substring(readFile1.length()-3,readFile1.length()).equals("bam")) {
    		isBAM = "true";
    	} else {
    		isBAM = "false";
    	}
    	pushFlowVariableString("readType", readType);
    	pushFlowVariableString("isBAM", isBAM);
    	
    	
		/**Create Output Table**/
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),}));
    	
    	DataCell[] c = new DataCell[]{
    			(FileCell) FileCellFactory.create(readFile1),
    			(FileCell) FileCellFactory.create(readFile2)};
    	
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
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String seq1 = m_readseqfile.getStringValue();
    	if(seq1.length() > 1) {
	    	if(!seq1.substring(seq1.length()-3,seq1.length()).equals("bam")) {
	    		if(!FileValidator.checkFastaFormat(seq1) && !FileValidator.checkFastqFormat(seq1)){
	                throw new InvalidSettingsException("Read sequences file no. 1 is not in Bam, FastQ or FastA format or does not contain nucleotide sequences!");
	        	}
	    	}
    	}
    	
    	if(m_readseqfile2.isEnabled()) {
    		String seq2 = m_readseqfile2.getStringValue();
    		if(seq2.length() > 1) {
	        	if(!FileValidator.checkFastaFormat(seq2) && !FileValidator.checkFastqFormat(seq2)){
	        		throw new InvalidSettingsException("Read sequences file no. 2 is not in FastQ or FastA format or does not contain nucleotide sequences!");
	        	}
    		}
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_readseqfile.saveSettingsTo(settings);
    	m_readseqfile2.saveSettingsTo(settings);
    	m_readType.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_readseqfile.loadSettingsFrom(settings);
    	m_readseqfile2.loadSettingsFrom(settings);
    	m_readType.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_readseqfile.validateSettings(settings);
    	m_readseqfile2.validateSettings(settings);
    	m_readType.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }

}

