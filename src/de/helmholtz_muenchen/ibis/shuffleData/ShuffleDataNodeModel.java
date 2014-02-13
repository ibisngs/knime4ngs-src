package de.helmholtz_muenchen.ibis.shuffleData;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Random;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
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
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

/**
 * This is the model implementation of ShuffleData.
 * Shuffle data within columns of the input data matrix.
 *
 * @author Jonas Zierer
 */
public class ShuffleDataNodeModel extends NodeModel {
	public static final String[] POSSIBLE_STRING_VALUES = {"String 1", "String 2", "Another String"};
	
	/** CFG KEYS */
	static final String CFGKEY_RANDOM_SEED    = "seed";
	static final String CFGKEY_COLUMNS        = "columns";
    
	/** SETTING MODELS */
    private final SettingsModelIntegerBounded m_seed  = new SettingsModelIntegerBounded(ShuffleDataNodeModel.CFGKEY_RANDOM_SEED, 5, Integer.MIN_VALUE, Integer.MAX_VALUE);
    private final SettingsModelFilterString   m_list   = new SettingsModelFilterString(ShuffleDataNodeModel.CFGKEY_COLUMNS);
    
    /**
     * Constructor for the node model.
     */
    protected ShuffleDataNodeModel() {
        super(1, 1);
    }

    private HashMap<String, Integer> getColIDs(DataTableSpec inData) throws InvalidSettingsException{
    	HashMap<String, Integer> colIDs = new HashMap<String, Integer>();
    	
    	for(String col: m_list.getIncludeList()){
    		int colID = inData.findColumnIndex(col);
    		
    		if(colID<0){
    			throw( new InvalidSettingsException("Column " + col + " does not exist in input data!")); 
    		}
    		colIDs.put(col, colID);
    	}
    	return(colIDs);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	HashMap<String, DataCell[]> colsToShuffle = new HashMap<String, DataCell[]>();
    	HashMap<String, Integer> ids = getColIDs(inData[0].getDataTableSpec());
    	String[] colNames = inData[0].getDataTableSpec().getColumnNames();
    	
    	// init data cell collection
    	for(String col: ids.keySet()){
    		colsToShuffle.put(col, new DataCell[inData[0].getRowCount()]);
    	}
    	
    	// collect all dataCells (which shall be shuffled
    	int rowCounter=0;
    	for(DataRow row: inData[0]){
    		for(String col: ids.keySet()){
    			colsToShuffle.get(col)[rowCounter] = row.getCell(ids.get(col));
    		}
    		rowCounter++;
    	}
    	
    	// shuffle data
    	Random random = new Random(m_seed.getIntValue());
    	int index;
    	DataCell temp;
    	for(String col: colsToShuffle.keySet()){
    		DataCell[] allCells = colsToShuffle.get(col);

    	    for (int i = allCells.length - 1; i > 0; i--){
    	        index = random.nextInt(i + 1);
    	        temp = allCells[index];
    	        allCells[index] = allCells[i];
    	        allCells[i] = temp;
    	    }
    		colsToShuffle.put(col, allCells);
    	}
    	
    	
    	// create output
    	rowCounter=0;
    	BufferedDataContainer output = exec.createDataContainer(inData[0].getDataTableSpec());
    	for(DataRow rowOld: inData[0]){
    		DataCell[] cellsNew = new DataCell[colNames.length];
    		for(int c=0; c<colNames.length; c++){
    			if(colsToShuffle.containsKey(colNames[c])){
    				cellsNew[c] = colsToShuffle.get(colNames[c])[rowCounter];
    			}else{
    				cellsNew[c] = rowOld.getCell(c);
    			}
    		}
    		
    		DataRow rowNew = new DefaultRow(rowOld.getKey(), cellsNew);
    		output.addRowToTable(rowNew);
    		rowCounter++;
    	}
    	output.close();
        return new BufferedDataTable[]{output.getTable()};
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
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return inSpecs;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_seed.saveSettingsTo(settings);
        m_list.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_seed.loadSettingsFrom(settings);
        m_list.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_seed.validateSettings(settings);
        m_list.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    }

}

