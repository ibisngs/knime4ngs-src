package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.modelExtraction;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsExtractionNodeModel extends RNodeModel {

	/** CFG KEYS */
	static final String CFGKEY_PERC_INCLUSION = "percentageInclusion";
	static final String CFGKEY_EV    = "e.v.";

    
	/** SETTING MODELS */
	private final SettingsModelDoubleBounded  m_percIncl = new SettingsModelDoubleBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_PERC_INCLUSION, 0.8, 0.0, 1.0);
    private final SettingsModelIntegerBounded m_ev       = new SettingsModelIntegerBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_EV, 5, 0, Integer.MAX_VALUE);

    
    /**
     * Constructor for the node model.
     */
    protected MixedGraphicalModelsExtractionNodeModel() {
        super(1, 2, "statistics" + File.separatorChar + "graphicalModels" + File.separatorChar + "mixedGraphicalModels.modelExtraction.R", new String[]{"--input"}, new String[]{"--edgeslist", "--adjacency"});
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
    	this.addArgument("--percIncl", m_percIncl.getDoubleValue());
    	this.addArgument("--ev"      , m_ev.getIntValue());

		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.getEdgeRankSpec(inData[0].getDataTableSpec())); // parse cell types
		out[1] = exec.createSpecReplacerTable(out[1], this.getAdjacencyMatSpec(inData[0].getDataTableSpec())); // parse cell types
		
		return(out);
	}
    
    
    private DataTableSpec getAdjacencyMatSpec(DataTableSpec inSpec) throws InvalidSettingsException{
    	ArrayList<String>   adjacencyColNames = new ArrayList<String>();
    	ArrayList<DataType> adjacencyColTypes = new ArrayList<DataType>();
    	HashMap<String, Boolean> nodeNames = new HashMap<String, Boolean>();
    	String[] inputColnames = inSpec.getColumnNames();
    	String[] colnamesSep;
    	for(String colname: inputColnames){
    		colnamesSep = colname.split(vertexSeparatorChar);
    		if(colnamesSep.length != 2){
    			throw new InvalidSettingsException("Can't figure out node names of edge >" +colname + "<\nNode names must not contain >"+vertexSeparatorChar+"<");
    		}
    		if(!nodeNames.containsKey(colnamesSep[0])){
    			nodeNames.put(colnamesSep[0], true);
    			adjacencyColNames.add(colnamesSep[0]);
    			adjacencyColTypes.add(DataType.getType(IntCell.class));
    		}
    		if(!nodeNames.containsKey(colnamesSep[1])){
    			nodeNames.put(colnamesSep[1], true);
    			adjacencyColNames.add(colnamesSep[1]);
    			adjacencyColTypes.add(DataType.getType(IntCell.class));
    		}
    	}
    	DataColumnSpec[] adjacencySpecs = DataTableSpec.createColumnSpecs(adjacencyColNames.toArray(new String[adjacencyColNames.size()]), adjacencyColTypes.toArray(new DataType[adjacencyColTypes.size()]));
    	
    	return(new DataTableSpec(adjacencySpecs));
	}
    
    private DataTableSpec getEdgeRankSpec(DataTableSpec inSpec){
    	DataColumnSpec[] edgeRankSpecs = DataTableSpec.createColumnSpecs(new String[]{"v1", "v2", "incl.n", "mean.rank", "incl.perc"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(StringCell.class), DataType.getType(IntCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(edgeRankSpecs));
	}
    
    
    public static final String vertexSeparatorChar = "~";
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
    	
    	
    	
    	
    	return new DataTableSpec[]{getEdgeRankSpec(inSpecs[0]), getAdjacencyMatSpec(inSpecs[0])};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        m_ev.saveSettingsTo(settings);
        m_percIncl.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
        m_ev.loadSettingsFrom(settings);
        m_percIncl.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_ev.validateSettings(settings);
        m_percIncl.validateSettings(settings);

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

