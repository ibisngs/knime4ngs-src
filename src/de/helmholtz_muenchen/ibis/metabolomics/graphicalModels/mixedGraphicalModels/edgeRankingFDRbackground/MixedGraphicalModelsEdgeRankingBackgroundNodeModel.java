package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRankingFDRbackground;

import java.io.File;
import java.io.IOException;
import java.util.Random;

import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.IntValue;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.port.PortObjectSpec;

import de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking.MixedGraphicalModelsEdgeRankingNodeFactory;
import de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking.MixedGraphicalModelsEdgeRankingNodeModel;
import de.helmholtz_muenchen.ibis.misc.shuffleData.ShuffleDataNodeModel;
import de.helmholtz_muenchen.ibis.utils.Global;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsEdgeRankingBackgroundNodeModel extends NodeModel {

	/** DEFAULT VALUES **/
	public static int DEFAULT_SAMPLENUM_FDR = 100;
	public static int INPORT_DATA                   = 0;
	public static int INPORT_RANDOM_SEEDS_STAB_SEL  = 1;
	public static int INPORT_RANDOM_SEEDS_EMP_PVALS = 2;
	
	/** CFG KEYS */
	public static final String CFGKEY_SS_SAMPLE_FDR           = "(empirical p-values) sample number";
	public static final String CFGKEY_RAND_SEED_EP_COL        = "(empirical p-values) random seed column";

	/** SETTING MODELS */
	private int m_sampleNum_FDR = 1;
	private String m_randomSeedsCol = MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS;
	
	private MixedGraphicalModelsEdgeRankingNodeModel edgeRanker;
	private MixedGraphicalModelsEdgeRankingNodeFactory edgeRankerFac;
	
	/**
	 * Constructor for the node model.
	 */
	protected MixedGraphicalModelsEdgeRankingBackgroundNodeModel() {
		super(Global.createOPOs(3,2,3), Global.createOPOs(2));
		edgeRankerFac = new MixedGraphicalModelsEdgeRankingNodeFactory();
		edgeRanker = edgeRankerFac.createNodeModel();
		
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws CanceledExecutionException{	
		////////////////////////////////////////////////////////////////////////////////////
		// RANDOM SEEDS FOR SHUFFLING
		////////////////////////////////////////////////////////////////////////////////////
		int randomSeeds[];		
		// GET RANDOM SEEDS FROM OPTIONAL INPUT
		BufferedDataTable randomSeedsTable = inData[INPORT_RANDOM_SEEDS_EMP_PVALS];
		if (randomSeedsTable != null){
			randomSeeds = new int[randomSeedsTable.getRowCount()];
			int randomSeedsColIdx = randomSeedsTable.getDataTableSpec().findColumnIndex(this.m_randomSeedsCol);
			if(randomSeedsColIdx == -1){
				throw new CanceledExecutionException("Can't find column >"+this.m_randomSeedsCol+"< in second input table!" ); 
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
			randomSeeds = new int[this.m_sampleNum_FDR];
			Random rand = new Random(System.currentTimeMillis());
			for(int i=0; i<randomSeeds.length; i++){
				randomSeeds[i] = rand.nextInt();
			}
		}
		exec.checkCanceled();

		////////////////////////////////////////////////////////////////////////////////////
		// DO EDGE RANKING
		////////////////////////////////////////////////////////////////////////////////////
		BufferedDataContainer edgeRankContainer = exec.createDataContainer(MixedGraphicalModelsEdgeRankingNodeModel.getEdgeRanksSpec(inData[0].getDataTableSpec()));
		BufferedDataContainer metaInfContainer = exec.createDataContainer(MixedGraphicalModelsEdgeRankingNodeModel.getMetaInfoSpec(inData[0].getDataTableSpec()));
		
		for(int i=0; i<randomSeeds.length; i++){
			// shuffle Data
			exec.checkCanceled();
			BufferedDataTable inputShuffled = ShuffleDataNodeModel.shuffleData(inData[INPORT_DATA], exec, inData[INPORT_DATA].getDataTableSpec().getColumnNames(), randomSeeds[i]);
			
			// edge Ranking
			BufferedDataTable[] result = this.edgeRanker.executeAnywhere(new BufferedDataTable[] {inputShuffled, inData[INPORT_RANDOM_SEEDS_STAB_SEL]}, exec);

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

			// TODO R- Logs
		}
		exec.checkCanceled();
		edgeRankContainer.close();
		metaInfContainer.close();

		return(new BufferedDataTable[]{edgeRankContainer.getTable(), metaInfContainer.getTable()});
	}

	
	protected PortObjectSpec[] configure(final PortObjectSpec[] inSpecs) throws InvalidSettingsException {
		PortObjectSpec[] datatableSpec = this.edgeRanker.configureAnywhere(new PortObjectSpec[]{inSpecs[INPORT_DATA], inSpecs[INPORT_RANDOM_SEEDS_STAB_SEL] });
		
		if(inSpecs[INPORT_RANDOM_SEEDS_EMP_PVALS] != null){
			DataTableSpec optionalInputSpec = (DataTableSpec)inSpecs[INPORT_RANDOM_SEEDS_EMP_PVALS];
			String[] validCols = Global.getValidCols(optionalInputSpec, IntValue.class);
			if(validCols.length==0){
				throw new InvalidSettingsException("No column with integer values in second input!");
			}
			if(this.m_randomSeedsCol.equals(MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS)){
				m_randomSeedsCol = validCols[0];
			}
		}else{
			m_randomSeedsCol = MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS;
		}
		return datatableSpec;
	}


	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		this.edgeRanker.resetAnywhere();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		this.edgeRanker.loadValidatedSettingsFromAnywhere(settings);
		
		// additional params
		this.m_sampleNum_FDR = settings.getInt(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_SS_SAMPLE_FDR, DEFAULT_SAMPLENUM_FDR);
		this.m_randomSeedsCol = settings.getString(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_RAND_SEED_EP_COL, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS);
		
	}

    /**
     * {@inheritDoc}
     */
	@Override
	 protected void saveSettingsTo(final NodeSettingsWO settings){
		this.edgeRanker.saveSettingsToAnywhere(settings);
		
		// additional params
		settings.addInt(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_SS_SAMPLE_FDR, this.m_sampleNum_FDR);
		settings.addString(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_RAND_SEED_EP_COL, this.m_randomSeedsCol);
	}

    /**
     * {@inheritDoc}
     */
	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {		
		this.edgeRanker.validateSettingsAnywhere(settings);
	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	this.edgeRanker.loadInternalsAnywhere(internDir, exec);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	this.edgeRanker.saveInternalsAnywhere(internDir, exec);
    }
}

