package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.edgeRanking;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.ParamUtils;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GraphicalModelsNodeModel extends RNodeModel {
	static final String[] GRAFO_RANKERS  = {"randomForest", "lasso", "randomForest.party", "regression"};
	static final String[] GRAFO_RANKTYPE = {"local", "global"};
	static final String NODE_ID_SEP = "~";
	
	static final String PARAMS_SEP     = ",";
	static final String PARAMS_SEP_2    = ";";
	static final String KEY_VALUE_SEP   = "=";
	static final String KEY_VALUE_SEP_2 = ":";

	static final int    DEFAULT_SS_SAMPLE_N = 100;
	static final double DEFAULT_SS_SAMPLE_S = 0.5;
	
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GraphicalModelsNodeModel.class);

	/** CFG KEYS */
	static final String CFGKEY_RANKER             = "rankers";
	static final String CFGKEY_PARAMS             = "rankerparams";
	static final String CFGKEY_VAR_TYPES          = "vartypes";
	static final String CFGKEY_SS_SAMPLE_N        = "stabselSamplenum";
	static final String CFGKEY_SS_RANKTYPE        = "stabselRankType";
	static final String CFGKEY_SS_SAMPLE_S        = "stabselSamplesize";
	static final String CFGKEY_RANDOM_SEED        = "seed";
	static final String CFGKEY_PARALLEL           = "parallel";

	/** SETTING MODELS */
	private final Map<String, String> m_ranker        = new HashMap<String, String>();
	private final HashMap<String, ArrayList<Pair<String,String>>> m_ranker_params = new HashMap<String, ArrayList<Pair<String,String>>>();
	private int m_stabSel_sampleNum = 1;
	private double m_stabSel_sampleSize = DEFAULT_SS_SAMPLE_S;
	private int m_rseed;
	private int m_parallel;
	
	private String m_rankerType = GRAFO_RANKTYPE[0];

	/**
	 * Constructor for the node model.
	 */
	protected GraphicalModelsNodeModel() {
		super(1, 2, "statistics" + File.separatorChar + "graphicalModels" + File.separatorChar + "runner.edgeranking.R", new String[]{"--data"}, new String[]{"--output", "--output2"});
		m_rseed = new Random().nextInt();
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
		this.addArgument("--classes"   , getVariableTypes(inData[0].getSpec()));

		// types
		Iterator<String> typesIt = getUniqueVariableTypes(inData[0].getSpec()).iterator();

		// Rankers + Params
		String ranker = "";
		String params = "";
		while(typesIt.hasNext()){
			String type = typesIt.next();
			ranker += type + KEY_VALUE_SEP + m_ranker.get(type) + PARAMS_SEP;
			params += type + KEY_VALUE_SEP_2 +ParamUtils.paramsToString(m_ranker_params.get(type), PARAMS_SEP, KEY_VALUE_SEP) + PARAMS_SEP_2; // TODO: use variables
		}
		this.addArgument("--ranker", ranker.substring(0, ranker.length()-  PARAMS_SEP.length()));
		this.addArgument("--param", params.substring( 0, params.length()-PARAMS_SEP_2.length()));
		
		// Other parameters

		this.addArgument("--sampleNum" , m_stabSel_sampleNum );
		this.addArgument("--sampleSize", m_stabSel_sampleSize);
		this.addArgument("--rankType"  , m_rankerType);
		this.addArgument("--rseed"     , m_rseed);
		this.addArgument("--cores"     , m_parallel);
		
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.edgeRanksSpec(inData[0].getDataTableSpec())); // pasre all to double cells
		
		return(out);
	}




	public static ArrayList<String> getVariableTypes(final DataTableSpec inSpecs){
		ArrayList<String> classes = new ArrayList<String>(inSpecs.getNumColumns());
		for(int c=0; c<inSpecs.getNumColumns(); c++){
			classes.add(c, inSpecs.getColumnSpec(c).getType().getPreferredValueClass().getName());
		}
		return(classes);
	}

	public static HashSet<String> getUniqueVariableTypes(final DataTableSpec inSpecs){
		HashSet<String> classes = new HashSet<String>();
		for(int c=0; c<inSpecs.getNumColumns(); c++){
			classes.add(inSpecs.getColumnSpec(c).getType().getPreferredValueClass().getName());
		}
		return(classes);
	}

	/**
	 * Passes the input spec to the output.
	 * 
	 * @param inSpecs The input spec.
	 * @return The generated output specs.
	 * @throws InvalidSettingsException If column to bin cannot be identified.
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		HashSet<String> types = getUniqueVariableTypes(inSpecs[0]);
		Iterator<String> typesIt = types.iterator();
		// add newly discovered types
		while(typesIt.hasNext()){
			String type = typesIt.next();
			if( ! this.m_ranker.containsKey(type)){
				this.m_ranker.put(type, GRAFO_RANKERS[0]);
				this.m_ranker_params.put(type, new ArrayList<Pair<String,String>>());
			}
		}
		return new DataTableSpec[]{edgeRanksSpec(inSpecs[0]),null};
	}

	private DataTableSpec edgeRanksSpec(DataTableSpec inSpecs){
		// create DataTableSpec for edge Ranks
		int nNodes = inSpecs.getNumColumns();
		String[] nodeNames = inSpecs.getColumnNames();
		int nEdges = nNodes*(nNodes-1)/2;
		// types
		DataType[] colTypes = new DataType[nEdges];
		String[] colNames   = new String[nEdges];
		for(int e = 0; e < nEdges; e++){
			colTypes[e] = DataType.getType(DoubleCell.class);
		}
		// edge names
		int eCounter=0;
		for(int j = 1; j < nNodes; j++) {
			for(int i = 0; i < j; i++) {
				colNames[eCounter++] = nodeNames[i] + NODE_ID_SEP + nodeNames[j];
			}
		}
		
		return(new DataTableSpec(colNames, colTypes));
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
		this.m_ranker.clear();
		this.m_ranker_params.clear();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		/* NORMAL SETTINGS */
		m_stabSel_sampleNum = settings.getInt(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_N, DEFAULT_SS_SAMPLE_N);
		m_stabSel_sampleSize = settings.getDouble(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_S, DEFAULT_SS_SAMPLE_S);
		m_rankerType = settings.getString(GraphicalModelsNodeModel.CFGKEY_SS_RANKTYPE , GraphicalModelsNodeModel.GRAFO_RANKTYPE[0]);
		m_rseed = settings.getInt(GraphicalModelsNodeModel.CFGKEY_RANDOM_SEED, new Random().nextInt());
		m_parallel = settings.getInt(GraphicalModelsNodeModel.CFGKEY_PARALLEL, Runtime.getRuntime().availableProcessors());
		
		/* ADVANCED TYPE/RANKER SETTINGS */
		String[] types       = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_VAR_TYPES   , new String[0]);
		String[] ranker      = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_RANKER      , new String[0]);
		String[] params      = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_PARAMS      , new String[0]);

		if(ranker.length != types.length){
			throw(new InvalidSettingsException("There has to be one ranker/ranker params entry per variable type!"));
		}

		for (int i = 0; i < types.length; i++) {
			m_ranker.put(types[i], ranker[i]);
			try {
				m_ranker_params.put(types[i], ParamUtils.stringToParams(params[i], PARAMS_SEP, KEY_VALUE_SEP));
			} catch (Exception e) {
				throw(new InvalidSettingsException("Can't parse parameters. " + e.getMessage()));
			}
		}
	}

	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings){
		/* NORMAL SETTINGS */
		settings.addInt(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_N, m_stabSel_sampleNum);
		settings.addDouble(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_S, m_stabSel_sampleSize);
		settings.addString(GraphicalModelsNodeModel.CFGKEY_SS_RANKTYPE, this.m_rankerType);
		settings.addInt(GraphicalModelsNodeModel.CFGKEY_RANDOM_SEED, m_rseed);
		settings.addInt(GraphicalModelsNodeModel.CFGKEY_PARALLEL, m_parallel);
		
		/* ADVANCED TYPE/RANKER SETTINGS */
		int nTypes = this.m_ranker.keySet().size();

		// init String arrays
		String[] types      = new String[nTypes];
		String[] ranker     = new String[nTypes];
		String[] params     = new String[nTypes];

		Iterator<String> typesIt = this.m_ranker.keySet().iterator();
		int tCounter=0;
		while(typesIt.hasNext()){
			String type = typesIt.next();
			types[tCounter]   = type;
			ranker[tCounter]  = this.m_ranker.get(type);
			params[tCounter]  = ParamUtils.paramsToString(this.m_ranker_params.get(type), PARAMS_SEP, KEY_VALUE_SEP);
			tCounter++;
		}
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_VAR_TYPES  , types);
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_RANKER     , ranker);
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_PARAMS     , params);
	}

	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {		
		/* NORMAL SETTINGS */
		if(settings.getInt(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_N, GraphicalModelsNodeModel.DEFAULT_SS_SAMPLE_N) < 1){
			throw(new InvalidSettingsException("Invalid number for samples for stability selection!")); 
		}
		
		if(settings.getDouble(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_S, GraphicalModelsNodeModel.DEFAULT_SS_SAMPLE_S)>1 | settings.getDouble(GraphicalModelsNodeModel.CFGKEY_SS_SAMPLE_S, GraphicalModelsNodeModel.DEFAULT_SS_SAMPLE_S)<0){
			throw(new InvalidSettingsException("Size of subsamples must be in [0.0;1.0]")); 
		}

		/* ADVANCED TYPE/RANKER SETTINGS */
		// TODO: validation of settings
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

