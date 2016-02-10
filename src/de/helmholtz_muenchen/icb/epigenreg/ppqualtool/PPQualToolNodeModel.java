package de.helmholtz_muenchen.icb.epigenreg.ppqualtool;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.ngs.*;


/**
 * This is the model implementation of PPQualTool.
 * 
 *
 * @author 
 */
public class PPQualToolNodeModel extends RNodeModel {
	
	//Static Variables:
	private static String PATH_PHANTOMPEAKQUALTOOL = "";
    
    // the logger instance
    private static final NodeLogger logger = NodeLogger
            .getLogger(PPQualToolNodeModel.class);
        
    /** the settings key which is used to retrieve and 
        store the settings (from the dialog or from a settings file)    
       (package visibility to be usable from the dialog). */
    
    static final String CFGKEY_OUTFILE = "Outfile";
	

    // example value: the models count variable filled from the dialog 
    // and used in the models execution method. The default components of the
    // dialog work with "SettingsModels".
    static final String DEFAULT_OUTFILE = "qc.txt";
    private final SettingsModelString SET_OUTFILE =
        new SettingsModelString(PPQualToolNodeModel.CFGKEY_OUTFILE, DEFAULT_OUTFILE);
    

    /**
     * Constructor for the node model.
     */
    protected PPQualToolNodeModel(){
    
        // TODO one incoming port and one outgoing port is assumed
    	//TODO add paths
        super(1, 1, PATH_PHANTOMPEAKQUALTOOL, new String[0], new String[0]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
        logger.info("Initializing R execution");
        
        //retrieve BAM file name from BufferedDataTabel
        String readsFile1 = inData[0].iterator().next().getCell(0).toString();
        
        //create command
        this.addArgument("-c", readsFile1);
        this.addArgument("-savp", "");
        this.addArgument("-out=", SET_OUTFILE.getStringValue());
        
        
        logger.info("");
        return(super.execute(inData, exec));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
        // TODO: check if user settings are available, fit to the incoming
        // table structure, and the incoming types are feasible for the node
        // to execute. If the node can execute in its current state return
        // the spec of its output data table(s) (if you can, otherwise an array
        // with null elements), or throw an exception with a useful user message
    	
    	//validate incoming data table spec
    	boolean hasBAMfile = false;
    	boolean bamCell = false;
    	//search filename (String) in input array
    	//we can only have one input DtatTableSpec and it has only one row and colum
    	//containing Bamfile name
    	//Normally, the is only one input incoming
    	 for (int i = 0; i < inSpecs[0].getNumColumns(); i++) {
             DataColumnSpec columnSpec = inSpecs[0].getColumnSpec(i);
             if(columnSpec.getType().isCompatible(BAMCell.class)){
            	 bamCell = true;
             }
             //TODO: Check if table name existst???
    	 }
    	 
    	 if(!bamCell) {
    		 throw new InvalidSettingsException("Input table must be a BAM file");
    	 }
    	 //Empty return, since no further processing planned
        return new DataTableSpec[]{new DataTableSpec()};
    }


 

    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // TODO load internal data. 
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
       
        // TODO save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }
    public void init(){
    	super.init();
    	addSetting(SET_OUTFILE);
    }

}

