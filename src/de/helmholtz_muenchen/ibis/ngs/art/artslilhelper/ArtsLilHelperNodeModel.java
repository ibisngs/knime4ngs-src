package de.helmholtz_muenchen.ibis.ngs.art.artslilhelper;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.StringUtils;
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

import de.helmholtz_muenchen.ibis.ngs.art.ArtNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;


/**
 * This is the model implementation of ArtsLilHelper.
 * 
 *
 * @author Tanzeem Haque
 */
public class ArtsLilHelperNodeModel extends NodeModel {
    
    // the logger instance
//    private static final NodeLogger logger = NodeLogger
//            .getLogger(ArtsLilHelperNodeModel.class);
//        
    /** the settings key which is used to retrieve and 
        store the settings (from the dialog or from a settings file)    
       (package visibility to be usable from the dialog). */
//	static final String CFGKEY_COUNT = "Count";
//
//    /** initial default count value. */
//    static final int DEFAULT_COUNT = 100;
//
//    // example value: the models count variable filled from the dialog 
//    // and used in the models execution method. The default components of the
//    // dialog work with "SettingsModels".
//    private final SettingsModelIntegerBounded m_count =
//        new SettingsModelIntegerBounded(ArtsLilHelperNodeModel.CFGKEY_COUNT,
//                    ArtsLilHelperNodeModel.DEFAULT_COUNT,
//                    Integer.MIN_VALUE, Integer.MAX_VALUE);
//    

	static boolean optionalPort = true;
	public static final String HG_19 = "/home/ibis/tanzeem.haque/Documents/hg19/"; 
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ArtNodeModel.class); //not used yet

    /**
     * Constructor for the node model.
     */
    protected ArtsLilHelperNodeModel() {
    
        // TODO one incoming port and one outgoing port is assumed
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	ArrayList<Chromosome_fastq> chr_fq_arrList = new ArrayList<>(24);
    	CloseableRowIterator it = inData[0].iterator();
//		String record = inData[0].getSpec().getName();
//		System.out.println(record);
    	String tmp_chr = "";
    	ArrayList<String> F_read = new ArrayList<>();
    	ArrayList<String> M_read = new ArrayList<>();
    	ArrayList<String> C_read = new ArrayList<>();
    	
    	int row_count = inData[0].getRowCount()-1;
		while (it.hasNext()) {
			DataRow row = it.next();
//			System.out.println(row.toString());
			String id = row.getKey().toString().substring("Row x. ".length());
			String chunk_idx = id.split("_")[0]; // 0_ 1_ 2_ ...24_
			String chr = id.split("_")[1]; // chr1, chr2 ... chr21
			String indiv = id.split("_")[2].substring(0, 1); // _F, _M, _C
			if (tmp_chr.equals("") && row_count == inData[0].getRowCount()-1)
				tmp_chr = chr;

			
			if (!tmp_chr.equals("")){

				/**
				 * The start
				 */
				if (!chr.equals(tmp_chr)){

					Chromosome_fastq chr_fq = new Chromosome_fastq(tmp_chr, 
							new Individual[]{new Individual("F", F_read),
							new Individual("M", M_read),
							new Individual("C", C_read)});
					chr_fq_arrList.add(chr_fq);

					F_read = new ArrayList<>();
					M_read = new ArrayList<>();
					C_read = new ArrayList<>();
				}
				if(indiv.equals("F")) { /*0 F_fwd F_rev*/
					F_read.add(chunk_idx+"\t"+row.getCell(0).toString()+"\t"+row.getCell(1).toString());
				}
				else if(indiv.equals("M")) {/*0 M_fwd M_rev*/
					M_read.add(chunk_idx+"\t"+row.getCell(0).toString()+"\t"+row.getCell(1).toString());
				}
				else if(indiv.equals("C")) {/*0 C_fwd C_rev*/
					C_read.add(chunk_idx+"\t"+row.getCell(0).toString()+"\t"+row.getCell(1).toString());
				}
				
				/**
				 * Extra for the last row
				 */
				if(row_count==0) {
					Chromosome_fastq chr_fq = new Chromosome_fastq(tmp_chr, 
							new Individual[]{new Individual("F", F_read),
							new Individual("M", M_read),
							new Individual("C", C_read)});
					chr_fq_arrList.add(chr_fq);
				}
			
			}

			if (row_count < inData[0].getRowCount()-1) 
				tmp_chr = chr;
			row_count--;


		}

		for (int i = 0; i < chr_fq_arrList.size(); i++) {
			String a = chr_fq_arrList.get(i).getChr() + "\t";
			
			for(int j = 0; j < chr_fq_arrList.get(i).getTrio_fastq().length; j++) { //getTrio_fastq is (Array of length 3) F->Arraylist of (chunk\tfwd\trev)
				String b = chr_fq_arrList.get(i).getTrio_fastq()[j].getIndiv()+"\t";
				
				mergeChunks(exec, 
						chr_fq_arrList.get(i).getTrio_fastq()[j].getFastq().size(), /*size of chunk*/
						chr_fq_arrList.get(i).getChr(), /*name of chr*/
						chr_fq_arrList.get(i).getTrio_fastq()[j].getIndiv()+"_fwd");/*individual and read type*/
				mergeChunks(exec, 
						chr_fq_arrList.get(i).getTrio_fastq()[j].getFastq().size(), /*size of chunk*/
						chr_fq_arrList.get(i).getChr(), /*name of chr*/
						chr_fq_arrList.get(i).getTrio_fastq()[j].getIndiv()+"_rev");
				for (int k = 0; k < chr_fq_arrList.get(i).getTrio_fastq()[j].getFastq().size(); k++) { //getFastq is an arraylist of strings
					String c = chr_fq_arrList.get(i).getTrio_fastq()[j].getFastq().get(k);						
					System.out.println(a + b + c);
				}
			}
		}
		
    	int col_num = 3; //ArtNodeModel.IDS.size();
    	DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];/**
    	 * create the column(s) for (multi)fastq output of a trio with 6 rows
    	 */
    	/**
    	for (int i = 0; i < col_num; i++) {
    	    allColSpecs[i] = new DataColumnSpecCreator("Col. " + (i+1) + ": " + OUT_COL + ArtNodeModel.IDS.get(i), FileCell.TYPE).createSpec();

    	}**/
    	allColSpecs[0] = new DataColumnSpecCreator("Col. 1: " + ArtNodeModel.OUT_COL + "fwd", FileCell.TYPE).createSpec();
    	allColSpecs[1] = new DataColumnSpecCreator("Col. 2: " + ArtNodeModel.OUT_COL + "rev", FileCell.TYPE).createSpec();
    	allColSpecs[2] = new DataColumnSpecCreator("Col. 3: " + "Reference fasta", FileCell.TYPE).createSpec();

    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
    	BufferedDataContainer container = exec.createDataContainer(outputSpec);
    	// let's add m_count rows to it
    	/**
    	 * creating the 6 rows
    	 * i = 0,1,2,3,4,5
    	 */
    	FileCell[] cells = new FileCell[col_num];

    	int row_idx = 0;
    	for (int i = 0; i < chr_fq_arrList.size(); i++) {
    		String chr = chr_fq_arrList.get(i).getChr();
        	pushFlowVariableString("Reference"+"_"+i, ArtsLilHelperNodeModel.HG_19 + chr + ".fa"); 

//    		System.out.println("IDs for output: " + chr);
    	    RowKey key = new RowKey("Row row row your boat");

    		for (int j = 0; j < 3; j++) {
    	        switch (j) {
    			case 0:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_F");
//    	            System.out.println("key 0: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_F_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_F_rev.fq");
    				cells[2] = (FileCell) FileCellFactory.create(ArtsLilHelperNodeModel.HG_19 + chr + ".fa");

    				DataRow row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			case 1:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_M");
//    	            System.out.println("key 1: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_M_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_M_rev.fq");
    				cells[2] = (FileCell) FileCellFactory.create(ArtsLilHelperNodeModel.HG_19 + chr + ".fa");

    				row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			case 2:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_C");
//    	            System.out.println("key 2: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_C_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_C_rev.fq");
    				cells[2] = (FileCell) FileCellFactory.create(ArtsLilHelperNodeModel.HG_19 + chr + ".fa");

    				row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			
    			}
//    	        row_idx++;

    		}    	    
    	    // check if the execution monitor was canceled
    	    exec.checkCanceled();   

        }
    	// once we are done, we close the container and return its table
        container.close();
    	BufferedDataTable out = container.getTable();
        return new BufferedDataTable[]{out};
    }

    private void mergeChunks(ExecutionContext exec, int size, String chr, String indiv_read) throws CanceledExecutionException, InterruptedException,
	ExecutionException, UnsuccessfulExecutionException, IOException {
		// TODO Auto-generated method stub
		ArrayList<String> merge_command = new ArrayList<String>(5);
		merge_command.add("sh /home/ibis/tanzeem.haque/Documents/Scripts/Sequenciator/mergeChunks.sh");
		merge_command.add((size-1) + "");
		merge_command.add(chr);
		merge_command.add(indiv_read);
		merge_command.add(ArtNodeModel.INTERNAL_OUTPUT_PATH);
		
		for (String s : merge_command) {
			System.out.print(s + " ");
		}
		System.out.println();
		
		/**Execute**/
		Executor.executeCommand(new String[]{StringUtils.join(merge_command, " ")},exec,LOGGER);
    	
    	for (int i = 0; i < size; i++) {
    		boolean fq_del = new File(ArtNodeModel.INTERNAL_OUTPUT_PATH + i+"_"+chr+"_"+indiv_read+".fq").delete();
        	System.out.println("FQDEL: " + fq_del);
    	}
    	
		
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
        
    	/**
    	 * Check OptionalInputPort
    	 */
		try{
			inSpecs[0].getColumnNames();
			String[] colNames = inSpecs[0].getColumnNames();
			
//			for (String s : colNames)

			if (colNames.length==2 && colNames[0].equals("Col. 1: " + ArtNodeModel.OUT_COL + "fwd") &&
					colNames[1].equals("Col. 2: " + ArtNodeModel.OUT_COL + "rev")) 
				optionalPort=true;
			else
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			

    	}catch(NullPointerException npe){
			optionalPort=false;
    	
    	}
		
    	int col_num = 3; //ArtNodeModel.IDS.size();
    	DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];
    	allColSpecs[0] = new DataColumnSpecCreator("Col. 1: " + ArtNodeModel.OUT_COL + "fwd", FileCell.TYPE).createSpec();
    	allColSpecs[1] = new DataColumnSpecCreator("Col. 2: " + ArtNodeModel.OUT_COL + "rev", FileCell.TYPE).createSpec();
    	allColSpecs[2] = new DataColumnSpecCreator("Col. 3: " + "Reference fasta", FileCell.TYPE).createSpec();

    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
    	String s = "Optional Port for Art's li'l helper is ";
    	s += (optionalPort)? "enabled": "diasbled";
		System.out.println(s);
        return new DataTableSpec[]{outputSpec};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        // TODO save user settings to the config object.
        
//        m_count.saveSettingsTo(settings);

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
        
//        m_count.loadSettingsFrom(settings);

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

//        m_count.validateSettings(settings);

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
    
    /**
     * 
     * 2 nested internal classes for the art inputs
     * the chromosome id along with the fastqs of 3 individuals
     * @author tanzeem.haque
     *
     */
    private static class Chromosome_fastq {
    	private String chr;
    	private Individual[] trio_fastq = new Individual[3];
    	
    	private Chromosome_fastq(String chr, Individual[] trio_fastq) {
    		setChr(chr);
    		setTrio_fastq(trio_fastq);
    	}

		public String getChr() {
			return chr;
		}

		public void setChr(String chr) {
			this.chr = chr;
		}

		public Individual[] getTrio_fastq() {
			return trio_fastq;
		}

		public void setTrio_fastq(Individual[] trio_fastq) {
			this.trio_fastq = trio_fastq;
		}
    }
    /**
     * 
     * A nested internal class in the nested internal class
     * @author tanzeem.haque
     *
     */
    private static class Individual {
    	private String indiv;
    	private ArrayList<String> fastq = new ArrayList<>();
    	
    	private Individual(String s, ArrayList<String> fastq) {
    		setIndiv(s);
    		setFastq(fastq);
    	}

		public String getIndiv() {
			return indiv;
		}

		public void setIndiv(String s) {
			this.indiv = s;
		}

		public ArrayList<String> getFastq() {
			return fastq;
		}

		public void setFastq(ArrayList<String> fastq) {
			this.fastq = fastq;
		}
	
    }

}

