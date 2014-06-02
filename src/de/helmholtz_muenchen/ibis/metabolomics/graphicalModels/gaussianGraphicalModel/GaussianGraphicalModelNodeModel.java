package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.gaussianGraphicalModel;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GaussianGraphicalModelNodeModel extends RNodeModel {

	/** CFG KEYS */
    
	/** SETTING MODELS */
 
    
    /**
     * Constructor for the node model.
     */
    protected GaussianGraphicalModelNodeModel() {
    	super(1, 1, "statistics" + File.separatorChar + "graphicalModels" + File.separatorChar + "gaussianGraphicalModels.R", new String[]{"--data"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{  	
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.getEdgeRankSpec(inData[0].getDataTableSpec())); // parse cell types

		return(out);
	}
    
    private DataTableSpec getEdgeRankSpec(DataTableSpec inSpec){
    	DataColumnSpec[] edgeRankSpecs = DataTableSpec.createColumnSpecs(new String[]{"v1", "v2", "pcor", "p", "q", "prob"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(StringCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(edgeRankSpecs));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	return new DataTableSpec[]{getEdgeRankSpec(inSpecs[0])};
    }

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
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec); // load output from stdout and stderr
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec); // save output from stdout and stderr
    }

}

