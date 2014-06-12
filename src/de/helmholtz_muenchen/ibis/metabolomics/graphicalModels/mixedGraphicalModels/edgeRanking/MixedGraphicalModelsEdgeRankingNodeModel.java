package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.IntValue;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.Global;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsEdgeRankingNodeModel extends RNodeModel {
	static final String[] GRAFO_RANKERS  = {"randomForest", "lasso", "randomForest.party", "regression"};
	static final String[] GRAFO_RANKTYPE = {"local", "global"};
	static final String NODE_ID_SEP = "~";
	
	static final String PARAMS_SEP     = ",";
	static final String PARAMS_SEP_2    = ";";
	static final String KEY_VALUE_SEP   = "=";
	static final String KEY_VALUE_SEP_2 = ":";

	static final int    DEFAULT_SS_SAMPLE_N = 100;
	static final double DEFAULT_SS_SAMPLE_S = 0.5;
	static final int    DEFAULT_CORES       = 1;
	
	protected static final NodeLogger LOGGER = NodeLogger.getLogger(MixedGraphicalModelsEdgeRankingNodeModel.class);

	/** CFG KEYS */
	static final String CFGKEY_RANKER             = "ranker";
	static final String CFGKEY_PARAMS             = "ranker parameters";
	static final String CFGKEY_VAR_TYPES          = "variable types";
	static final String CFGKEY_PARALLEL           = "parallel cores";
	static final String CFGKEY_SS_SAMPLE_N        = "(stability selection) sample number";
	static final String CFGKEY_SS_SAMPLE_S        = "(stability selection) sample size";
	static final String CFGKEY_SS_RANKTYPE        = "(stability selection) ranking type";
	static final String CFGKEY_RANDOM_SEED_COL    = "(stability selection) random seed column";

	/** SETTING MODELS */
	private final Map<String, String> m_ranker        = new HashMap<String, String>();
	private final HashMap<String, ArrayList<Pair<String,String>>> m_ranker_params = new HashMap<String, ArrayList<Pair<String,String>>>();
	private int m_stabSel_sampleNum     = DEFAULT_SS_SAMPLE_N;
	private double m_stabSel_sampleSize = DEFAULT_SS_SAMPLE_S;
	private String m_rseed_col;
	private int m_cores              = DEFAULT_CORES;
	private String m_rankerType         = GRAFO_RANKTYPE[0];
	public static String DEFAULT_AVAILABLE_COLS = "2nd input not available";

	/**
	 * Constructor for the node model.
	 */
	protected MixedGraphicalModelsEdgeRankingNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes ) {
		super(inPortTypes, outPortTypes, "statistics" + File.separatorChar + "graphicalModels" + File.separatorChar + "mixedGraphicalModels.edgeranking.R", new String[]{"--data"}, new String[]{"--output", "--output2"});
	}
	
	
	 /**
	  * Constructor for the node model.
	  */
	protected MixedGraphicalModelsEdgeRankingNodeModel() {
		this(Global.createOPOs(2,2), Global.createOPOs(2));
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		////////////////////////////////////////////////////////////////////////////////////
		// WRITE VARIABLE CLASSES
		////////////////////////////////////////////////////////////////////////////////////
		File tmpFile = null;
		try {
			tmpFile = File.createTempFile("knime_R_connector_" + this.SCRIPT.replaceAll(File.separatorChar+"", "_") + "_input_" , ".csv");
			ArrayList<String> classes = Global.getVariableTypes(inData[0].getDataTableSpec());
			StringBuffer classesFileContent = new StringBuffer();
			classesFileContent.append("\"" + StringUtils.join(inData[0].getDataTableSpec().getColumnNames(), "\";\"") + "\"");
			classesFileContent.append("\n");
			classesFileContent.append("\"" + StringUtils.join(classes, "\";\"") + "\"");
			classesFileContent.append("\n");
			FileUtils.write(tmpFile, classesFileContent.toString());
			this.addArgument("--fclasses", tmpFile.getCanonicalPath());
		} catch (IOException e) {
			LOGGER.error("Unable to create temp file!");
			throw(new IOException("unable to create temp file!" + e.getMessage()));
		}

		
		////////////////////////////////////////////////////////////////////////////////////
		// RANDOM SEEDS FOR SAMPLING
		////////////////////////////////////////////////////////////////////////////////////
		int randomSeeds[];		
		// GET RANDOM SEEDS FROM OPTIONAL INPUT
		BufferedDataTable randomSeedsTable = inData[1];
		if (randomSeedsTable != null){
			randomSeeds = new int[randomSeedsTable.getRowCount()];
			int randomSeedsColIdx = randomSeedsTable.getDataTableSpec().findColumnIndex(this.m_rseed_col);
			if(randomSeedsColIdx == -1){
				throw new InvalidSettingsException("Can't find column >"+this.m_rseed_col+"< in second input table!" ); 
			}
			CloseableRowIterator rit = randomSeedsTable.iterator();
			int rowCounter=0;
			while(rit.hasNext()){
				DataRow row = rit.next();
				IntValue cell = (IntValue)row.getCell(randomSeedsColIdx);
				randomSeeds[rowCounter++] = cell.getIntValue();
			}
		// GENERATE RANDOM SEEDS
		}else{
			randomSeeds = new int[this.m_stabSel_sampleNum];
			Random rand = new Random(System.currentTimeMillis());
			for(int i=0; i<randomSeeds.length; i++){
				randomSeeds[i] = rand.nextInt();
			}
		}
		
		
		////////////////////////////////////////////////////////////////////////////////////
		// PARSE PARAMETERS
		////////////////////////////////////////////////////////////////////////////////////		
		// types
		Iterator<String> typesIt = Global.getUniqueVariableTypes(inData[0].getSpec()).iterator();

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
		this.addArgument("--sampleSize", m_stabSel_sampleSize);
		this.addArgument("--ranktype"  , m_rankerType);
		this.addArgument("--cores"     , m_cores);
		
		
		////////////////////////////////////////////////////////////////////////////////////
		// DO EDGE RANKING
		////////////////////////////////////////////////////////////////////////////////////
		BufferedDataContainer edgeRankContainer = exec.createDataContainer(getEdgeRanksSpec(inData[0].getDataTableSpec()));
		BufferedDataContainer metaInfContainer = exec.createDataContainer(getMetaInfoSpec(inData[0].getDataTableSpec()));
		
		
		for(int i=0; i<randomSeeds.length; i++){
			exec.checkCanceled();
			this.addArgument("--rseed", randomSeeds[i]);
			this.addArgument("--sampleNum" , 1 );
			BufferedDataTable[] result = super.execute(new BufferedDataTable[]{inData[0]}, exec);
			// rename edgeRank rows
			exec.checkCanceled();
			CloseableRowIterator it = result[0].iterator();
			while(it.hasNext()){
				DataRow row = it.next();
				edgeRankContainer.addRowToTable(new DefaultRow(row.getKey().getString() + "."+i, row));
			}
			
			// rename metainfo rows
			exec.checkCanceled();
			it = result[1].iterator();
			while(it.hasNext()){
				DataRow row = it.next();
				metaInfContainer.addRowToTable(new DefaultRow(row.getKey().getString() + "."+i, row));
			}
		}
		exec.checkCanceled();
		edgeRankContainer.close();
		metaInfContainer.close();
		return(new BufferedDataTable[]{edgeRankContainer.getTable(), metaInfContainer.getTable()});
	}

	/**
	 * Passes the input spec to the output.
	 * 
	 * @param inSpecs The input spec.
	 * @return The generated output specs.
	 * @throws InvalidSettingsException If column to bin cannot be identified.
	 */
	@Override
	protected PortObjectSpec[] configure(final PortObjectSpec[] inSpecs) throws InvalidSettingsException {
		HashSet<String> types = Global.getUniqueVariableTypes((DataTableSpec)inSpecs[0]);
		Iterator<String> typesIt = types.iterator();
		// add newly discovered types
		while(typesIt.hasNext()){
			String type = typesIt.next();
			if( ! this.m_ranker.containsKey(type)){
				this.m_ranker.put(type, GRAFO_RANKERS[0]);
				this.m_ranker_params.put(type, new ArrayList<Pair<String,String>>());
			}
		}
		
		return new PortObjectSpec[]{getEdgeRanksSpec((DataTableSpec)inSpecs[0]),getMetaInfoSpec((DataTableSpec)inSpecs[0])};
	}
	
	
	public static DataTableSpec getMetaInfoSpec(DataTableSpec inSpec){
		DataColumnSpec[] metaInfoSpecs = DataTableSpec.createColumnSpecs(new String[]{"variables", "stabSel.Sample.num", "n", "p", "seed", "time"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(IntCell.class), DataType.getType(IntCell.class), DataType.getType(IntCell.class), DataType.getType(IntCell.class), DataType.getType(DoubleCell.class)});
		return(new DataTableSpec(metaInfoSpecs));
	}
	
	public static DataTableSpec getEdgeRanksSpec(DataTableSpec inSpec){
		// create DataTableSpec for edge Ranks
		int nNodes = inSpec.getNumColumns();
		String[] nodeNames = inSpec.getColumnNames();
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
		m_stabSel_sampleNum = settings.getInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_N, DEFAULT_SS_SAMPLE_N);
		m_stabSel_sampleSize = settings.getDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, DEFAULT_SS_SAMPLE_S);
		m_rankerType = settings.getString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_RANKTYPE , MixedGraphicalModelsEdgeRankingNodeModel.GRAFO_RANKTYPE[0]);
		m_rseed_col = settings.getString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANDOM_SEED_COL, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS);
		m_cores     = settings.getInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARALLEL, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_CORES);
		/* ADVANCED TYPE/RANKER SETTINGS */
		String[] types       = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_VAR_TYPES   , new String[0]);
		String[] ranker      = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANKER      , new String[0]);
		String[] params      = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARAMS      , new String[0]);

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
		settings.addInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_N, m_stabSel_sampleNum);
		settings.addDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, m_stabSel_sampleSize);
		settings.addString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_RANKTYPE, this.m_rankerType);
		settings.addString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANDOM_SEED_COL, m_rseed_col);
		settings.addInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARALLEL, m_cores);
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
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_VAR_TYPES  , types);
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANKER     , ranker);
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARAMS     , params);
	}

	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {		
		/* NORMAL SETTINGS */
		if(settings.getInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_N, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_N) < 1){
			throw(new InvalidSettingsException("Invalid number for samples for stability selection!")); 
		}
		
		if(settings.getDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_S)>1 | settings.getDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_S)<0){
			throw(new InvalidSettingsException("Size of subsamples must be in [0.0;1.0]")); 
		}

		/* ADVANCED TYPE/RANKER SETTINGS */
		// TODO: validation of settings
	}

	public BufferedDataTable[] executeAnywhere(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		return(this.execute(inData, exec));
	}
	public void resetAnywhere() {
		this.reset();
	}
    public void validateSettingsAnywhere(NodeSettingsRO settings) throws InvalidSettingsException {		
    	this.validateSettings(settings);
    }
    public void loadValidatedSettingsFromAnywhere(NodeSettingsRO settings) throws InvalidSettingsException {		
    	this.loadValidatedSettingsFrom(settings);
    }
    public void saveSettingsToAnywhere(NodeSettingsWO settings) {		
    	this.saveSettingsTo(settings);
    }
    public PortObjectSpec[] configureAnywhere(final PortObjectSpec[] inSpecs) throws InvalidSettingsException {
    	return(this.configure(inSpecs));
    }
    public void loadInternalsAnywhere(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	this.loadInternals(internDir, exec);
    }
    public void saveInternalsAnywhere(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	this.saveInternals(internDir, exec);
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

