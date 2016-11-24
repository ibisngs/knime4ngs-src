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


package de.helmholtz_muenchen.ibis.ngs.vcfloader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * This is the model implementation of VCFLoader.
 * 
 *
 * @author Marie-Sophie Friedl
 */
public class VCFLoaderNodeModel extends NodeModel {
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger logger = NodeLogger.getLogger(VCFLoaderNodeModel.class);
    
    //path to first vcf file
    static final String CFGKEY_VCF1 ="vcf1";
    static final String DEF_VCF1="";
    private final SettingsModelString m_vcf1 = new SettingsModelString(CFGKEY_VCF1, DEF_VCF1);
    // type of variants in first vcf file
    static final String CFGKEY_TYPE1="type1";
    static final String[] AVAIL_TYPES={"SNPs", "Indels", "Insertions", "Deletions"};
    static final String DEF_TYPE1=AVAIL_TYPES[0];
    private final SettingsModelString m_type1 = new SettingsModelString(CFGKEY_TYPE1, DEF_TYPE1);
    // upload second file
    static final String CFGKEY_SECONDFILE="secondfile";
    static final boolean DEF_SECONDFILE=false;
    private final SettingsModelBoolean m_secondfile = new SettingsModelBoolean(CFGKEY_SECONDFILE, DEF_SECONDFILE);
    // path to second vcf file
    static final String CFGKEY_VCF2="vcf2";
    static final String DEF_VCF2="";
    private final SettingsModelString m_vcf2 = new SettingsModelString(CFGKEY_VCF2, DEF_VCF2);
    // type of variants in second vcf file
    static final String CFGKEY_TYPE2="type2";
    static final String DEF_TYPE2=AVAIL_TYPES[1];
    private final SettingsModelString m_type2 = new SettingsModelString(CFGKEY_TYPE2, DEF_TYPE2);

    

    /**
     * Constructor for the node model.
     */
    protected VCFLoaderNodeModel() {
    
        // no inport, one outport
        super(0, 1);
        
        m_type2.setEnabled(false);
        m_vcf2.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//check node configurations
    	
    	//check first vcf file
        if(m_vcf1.getStringValue().equals("")){
        	throw new Exception("Missing path to vcf file: You have to configure the node before executing it!");
        }
        if(!Files.exists(Paths.get(m_vcf1.getStringValue()))){
        	throw new Exception("Path to vcf file: "+m_vcf1.getStringValue()+" does not exist!");
        }
        
        //check second file
        if(m_secondfile.getBooleanValue()){
        	
        	if(m_type1.getStringValue().equals(m_type2.getStringValue())){
        		throw new Exception("Variant type of both vcf files has to be different!");
        	}
        	
            if(m_vcf2.getStringValue().equals("")){
            	throw new Exception("Missing path to second vcf file: You have to configure the node before executing it!");
            }
            if(!Files.exists(Paths.get(m_vcf2.getStringValue()))){
            	throw new Exception("Path to second vcf file: "+m_vcf2.getStringValue()+" does not exist!");
            }        	
        }

        //create output table
        
    	// number of output columns
    	int colcount=1;
    	if(m_secondfile.getBooleanValue()){
        	colcount=2;
    	}
    	
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[colcount];
    	
		//SNPs
    	if(m_type1.getStringValue().equals(AVAIL_TYPES[0])){
    		colspec[0]=new DataColumnSpecCreator("Path2VCFsnpFile", StringCell.TYPE).createSpec();
    	}
    	//Indels
    	else if(m_type1.getStringValue().equals(AVAIL_TYPES[1])){
    		colspec[0]=new DataColumnSpecCreator("Path2VCFindelFile", StringCell.TYPE).createSpec();
    	}
    	//insertions
    	else if(m_type1.getStringValue().equals(AVAIL_TYPES[2])){
    		colspec[0]=new DataColumnSpecCreator("Path2VCFinsertionsFile", StringCell.TYPE).createSpec();
    	}
    	//deletions
    	else if(m_type1.getStringValue().equals(AVAIL_TYPES[3])){
    		colspec[0]=new DataColumnSpecCreator("Path2VCFdeletionsFile", StringCell.TYPE).createSpec();
    	}
	    	
    	if(m_secondfile.getBooleanValue()){
    		//SNPs
	    	if(m_type2.getStringValue().equals(AVAIL_TYPES[0])){
	    		colspec[1]=new DataColumnSpecCreator("Path2VCFsnpFile", StringCell.TYPE).createSpec();
	    	}
	    	//Indels
	    	else if(m_type2.getStringValue().equals(AVAIL_TYPES[1])){
	    		colspec[1]=new DataColumnSpecCreator("Path2VCFindelFile", StringCell.TYPE).createSpec();
	    	}
	    	//insertions
	    	else if(m_type2.getStringValue().equals(AVAIL_TYPES[2])){
	    		colspec[1]=new DataColumnSpecCreator("Path2VCFinsertionsFile", StringCell.TYPE).createSpec();
	    	}
	    	//deletions
	    	else if(m_type2.getStringValue().equals(AVAIL_TYPES[3])){
	    		colspec[1]=new DataColumnSpecCreator("Path2VCFdeletionsFile", StringCell.TYPE).createSpec();
	    	}
    	}

  	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    
	    // fill string cells
	    DataCell [] row = new DataCell [colcount];
	    row[0]=new StringCell(m_vcf1.getStringValue());
	    if(m_secondfile.getBooleanValue()){
	    	row[1]=new StringCell(m_vcf2.getStringValue());
	    }

	    
	    //create row and add it to the container
	    c.addRowToTable(new DefaultRow("row0", row));
	    
	    //create final table
	    c.close();
	    BufferedDataTable out=c.getTable();
 
        return new BufferedDataTable [] {out};

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
        // check if user settings are available, fit to the incoming
        // table structure, and the incoming types are feasible for the node
        // to execute. If the node can execute in its current state return
        // the spec of its output data table(s) (if you can, otherwise an array
        // with null elements), or throw an exception with a useful user message

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        m_vcf1.saveSettingsTo(settings);
        m_type1.saveSettingsTo(settings);
        m_secondfile.saveSettingsTo(settings);
        m_vcf2.saveSettingsTo(settings);
        m_type2.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	m_vcf1.loadSettingsFrom(settings);
    	m_type1.loadSettingsFrom(settings);
    	m_secondfile.loadSettingsFrom(settings);
    	m_vcf2.loadSettingsFrom(settings);
    	m_type2.loadSettingsFrom(settings);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	m_vcf1.validateSettings(settings);
    	m_type1.validateSettings(settings);
    	m_secondfile.validateSettings(settings);
    	m_vcf2.validateSettings(settings);
    	m_type2.validateSettings(settings);
    	
    	//check if files have same variant type
    	if(settings.getBoolean(CFGKEY_SECONDFILE) && settings.getString(CFGKEY_TYPE1).equals(settings.getString(CFGKEY_TYPE2))){
    		throw new InvalidSettingsException("Variant type of both vcf files has to be different!");
    	}

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }

}

