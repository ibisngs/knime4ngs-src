package de.helmholtz_muenchen.ibis.ngs.frostslilhelper;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.container.CloseableRowIterator;
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

import de.helmholtz_muenchen.ibis.ngs.frost.FrostNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.frost.FrostRunner;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;


/**
 * This is the model implementation of Frostslilhepler.
 * 
 *
 * @author 
 */
public class FrostsLilHelperNodeModel extends NodeModel {
    
	public final static int chunk_length = 10000000;
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FrostsLilHelperNodeModel.class); //not used yet

    

    /**
     * Constructor for the node model.
     */
    protected FrostsLilHelperNodeModel() {
    
        // TODO one incoming port and one outgoing port is assumed
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	CloseableRowIterator it = inData[0].iterator();
    	int row_count = inData[0].getRowCount()-1;
		while (it.hasNext()) {
			DataRow row = it.next();
			String id = row.getKey().toString().substring("Row: ".length());
			String[] trio = {
					row.getCell(0).toString(), //F_0
					row.getCell(1).toString(), //F_1
					row.getCell(2).toString(), //M_0
					row.getCell(3).toString(), //M_1
					row.getCell(4).toString(), //C_0
					row.getCell(5).toString() // C_1
					};
			
			for (int i = 0; i < trio.length; i++) {
				FrostRunner.getChunks(trio[i]);
			}
		}

		return output(exec);
        
    }

    private BufferedDataTable[] output(final ExecutionContext exec)
			throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException, IOException {
		
    	int col_num = 6; // 6 fastas
        DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];
        
        /**
         * create the column(s) for (multi)fasta output of a trio
         */ 
        allColSpecs[0] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M0.fa", FileCell.TYPE).createSpec();
		allColSpecs[1] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M1.fa", FileCell.TYPE).createSpec();
		allColSpecs[2] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F0.fa", FileCell.TYPE).createSpec();
		allColSpecs[3] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F1.fa", FileCell.TYPE).createSpec();
		allColSpecs[4] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C0.fa", FileCell.TYPE).createSpec();
		allColSpecs[5] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C1.fa", FileCell.TYPE).createSpec();
			
		
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
        BufferedDataContainer container = exec.createDataContainer(outputSpec);       
        /**
         * creating the rows
         * i = 0,1, ... #IDS
         */  
		FileCell[] cells = new FileCell[6];
//		System.out.println("ID LIST now: " + FrostRunner.id_list.size());
		for (int j = 0; j < FrostRunner.id_list.size(); j++) {
			RowKey key = new RowKey("Row: " + FrostRunner.id_list.get(j));
//			System.out.println(key.toString());
			
			for (int i = 0; i < 6 ; i++) {
			// the cells of the current row, the types of the cells must match
			// number of cells = col_num
			// the column spec (see above)

				switch (i) {
					case 0:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_M0.fa");
						break;
					case 1:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_M1.fa");
						break;
					case 2:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_F0.fa");
						break;
					case 3:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_F1.fa");
						break;
					case 4:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_C0.fa");
						break;
					case 5:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_C1.fa");
						break;
				}
					
			}
			DataRow row = new DefaultRow(key, cells);
            container.addRowToTable(row);
		}
		FrostRunner.id_list.clear();
		// check if the execution monitor was canceled
        exec.checkCanceled();
		// once we are done, we close the container and return its table
        container.close();
        BufferedDataTable out = container.getTable();
        return new BufferedDataTable[]{out};
        

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

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        // TODO save user settings to the config object.
        

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO load (valid) settings from the config object.
        // It can be safely assumed that the settings are valided by the 
        // method below.
        

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO check if the settings could be applied to our model
        // e.g. if the count is in a certain range (which is ensured by the
        // SettingsModel).
        // Do not actually set any values of any member variables.


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

}

