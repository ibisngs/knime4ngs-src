package de.helmholtz_muenchen.ibis.CountMissing;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;


/**
 * This is the model implementation of CountMissingNodeModel.
 *
 * @author Jonas Zierer
 */
public class CountMissingNodeModel extends NodeModel {

	/** CFG KEYS */
	static final String CFGKEY_PERCENT   = "de.helmholtz_muenchen.ibis.CountMissing.percent";
    
	/** SETTING MODELS */
    private final SettingsModelBoolean m_percent   = new SettingsModelBoolean(CountMissingNodeModel.CFGKEY_PERCENT, true);

    
    /**
     * Constructor for the node model.
     */
    protected CountMissingNodeModel() {
        super(1, 2);       
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws CanceledExecutionException{
    	BufferedDataTable data = inData[0];
    	
    	// init objects
    	int rowsTotal  = data.getRowCount();
    	int colsTotal  = inData[0].getDataTableSpec().getNumColumns();
    	int rowCounter = 0;
    	int colCounter = 0;
    	String[] rowNames = new String[rowsTotal];
    	CloseableRowIterator it_row = data.iterator();
    	DataRow r;
    	DataCell cell;
    	
    	// create output
    	int[] numberOfMissingsCol = new int[colsTotal];
    	int[] numberOfMissingsRow = new int[rowsTotal];
    	
    	// check if missing for each cell
    	while(it_row.hasNext()){
    		// test if still executing
    		try{
    			exec.checkCanceled();
    		} catch(CanceledExecutionException e){
    			it_row.close();
    			throw(e);
    		}
    		exec.setProgress((double)rowCounter/rowsTotal);
    		
    		// process row
    		colCounter = 0;
    		r = it_row.next();
    		rowNames[rowCounter] = r.getKey().getString();
    		Iterator<DataCell> it_cell = r.iterator();
    		while(it_cell.hasNext()){
    			cell = it_cell.next();
    			if(cell.isMissing()){
					numberOfMissingsCol[colCounter]++;
					numberOfMissingsRow[rowCounter]++;
    			}
    			colCounter++;
    		}
    		rowCounter++;
    	}

   
    	// convert result
    	return(new BufferedDataTable[]{missingsCol(inData[0], exec, numberOfMissingsCol), missingsRow(inData[0], exec, numberOfMissingsRow, rowNames)});	

	}
    
    
	//////////////////////////////////////////////////////////////////////////
	// PARSE ARRAYS TO DATATABLES
	//////////////////////////////////////////////////////////////////////////
    private BufferedDataTable missingsRow(final BufferedDataTable inData, final ExecutionContext exec, int[] missingValues, String[] rownames){
    	BufferedDataContainer container = exec.createDataContainer(getRowTableSpec(inData.getDataTableSpec()));
    	for(int i=0; i<missingValues.length; i++){
            RowKey key = new RowKey(rownames[i]); // TODO row key from input table
            DataCell[] cells = new DataCell[1];
            if(m_percent.getBooleanValue()){
            	cells[0] = new DoubleCell((double)missingValues[i]/inData.getDataTableSpec().getNumColumns());
            }else{
            	cells[0] = new IntCell(missingValues[i]);
            }
            DataRow row = new DefaultRow(key, cells);
            container.addRowToTable(row);
    	}
    	container.close();
    	return container.getTable();
    }
    
    private BufferedDataTable missingsCol(final BufferedDataTable inData, final ExecutionContext exec, int[] missingValues){
    	BufferedDataContainer container = exec.createDataContainer(getColTableSpec(inData.getDataTableSpec()));
    	RowKey key = new RowKey("Missing Values"); // TODO row key from input table
    	DataCell[] cells = new DataCell[missingValues.length];
    	for(int c=0; c < missingValues.length; c++){
    		if(m_percent.getBooleanValue()){
    			cells[c] = new DoubleCell((double)missingValues[c]/inData.getRowCount());
    		}else{
    			cells[c] = new IntCell(missingValues[c]);	
    		}
    	}
    	DataRow row = new DefaultRow(key, cells);
    	container.addRowToTable(row);
    	container.close();
    	return(container.getTable());
    }
    
    
    
	//////////////////////////////////////////////////////////////////////////
	// CREATE OUTPUT TABLESPECS
	//////////////////////////////////////////////////////////////////////////
    private DataTableSpec getColTableSpec(final DataTableSpec inSpecs){
    	String[] colnames = inSpecs.getColumnNames();
    	int colNum        = inSpecs.getNumColumns();
        DataColumnSpec[] allColSpecs = new DataColumnSpec[colNum];
        
        for(int c=0; c<colNum; c++){
        	if(m_percent.getBooleanValue()){
        		allColSpecs[c] = new DataColumnSpecCreator(colnames[c], DoubleCell.TYPE).createSpec();
        	}else{
        		allColSpecs[c] = new DataColumnSpecCreator(colnames[c], IntCell.TYPE).createSpec();
        	}
        }
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
        
        return(outputSpec);
    }
    
    private DataTableSpec getRowTableSpec(final DataTableSpec inSpecs){
    	DataColumnSpec[] allColSpecs = new DataColumnSpec[1];
    	if(m_percent.getBooleanValue()){
    		allColSpecs[0] = new DataColumnSpecCreator("Missing Values", DoubleCell.TYPE).createSpec();
    	}else{
    		allColSpecs[0] = new DataColumnSpecCreator("Missing Values", IntCell.TYPE).createSpec();
    	}
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
        return(outputSpec);
    }
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return new DataTableSpec[]{getColTableSpec(inSpecs[0]), getRowTableSpec(inSpecs[0])};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_percent.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_percent.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_percent.validateSettings(settings);
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

	@Override
	protected void reset() {		
	}

}

