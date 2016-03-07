package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.GVCFCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

/**
 * This is the model implementation of FileLoader.
 * 
 *
 * @author Tim Jeske
 */
public class FileLoaderNodeModel extends NodeModel {
	
	public static final String CFGKEY_INFILE1 = "infile1";
	public static final String CFGKEY_INFILE2 = "infile2";
	
	private final SettingsModelString m_infile1 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE1,"");
	private final SettingsModelString m_infile2 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE2,"");
	
	public static final String OUT_COL1 = "Path2File1";
	public static final String OUT_COL2 = "Path2File2";
	
	private static final String [] ENDINGS = {"",".vcf",".gvcf",".fastq",".fq",".bam",".sam"};
	private static final DataType [] TYPES = {FileCell.TYPE, VCFCell.TYPE, GVCFCell.TYPE, FastQCell.TYPE, FastQCell.TYPE, BAMCell.TYPE, SAMCell.TYPE};
	
	private boolean secondOk = false;
	private boolean is_list = false;
	private String sep;
	
	private DataColumnSpec dcs1 = null;
	private DataColumnSpec dcs2 = null;
	private DataColumnSpec [] specs = null;
	
    /**
     * Constructor for the node model.
     */
    protected FileLoaderNodeModel() {
    	super(0, 1);
    	m_infile2.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	ArrayList<String> in1_list = new ArrayList<>();
    	ArrayList<String> in2_list = new ArrayList<>();
    	
    	String in1 = m_infile1.getStringValue();
    	String in2 = m_infile2.getStringValue();
    	
    	if(is_list) {
    		BufferedReader bw = Files.newBufferedReader(Paths.get(in1));
    		String line;
    		while((line=bw.readLine())!=null) {
    			if(line.trim().equals("")) continue;
    			if(secondOk) {
    				String [] col = line.split(this.sep);
    				in1_list.add(col[0]);
    				if(col.length<2) {
    					setWarningMessage("Second column contains less entries than first! Ignoring those lines!");
    				}
    				in2_list.add(col[1]);
    			} else {
    				in1_list.add(line.trim());
    			}
    		}
    		
    		bw.close();
    	} else {
    		in1_list.add(in1);
    		in2_list.add(in2);
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(specs));
    	
    	for(int i = 0; i < in1_list.size(); i++) {
    		FileCell [] fileCell = new FileCell[] {
        			(FileCell) FileCellFactory.create(in1_list.get(i))};
        	
        	if(secondOk) {
        		fileCell = new FileCell[] {
            			(FileCell) FileCellFactory.create(in1_list.get(i)),
            			(FileCell) FileCellFactory.create(in2_list.get(i))};
        	}
        	FileCell[] c = fileCell;
    		cont.addRowToTable(new DefaultRow("Row"+i,c));
    	}
    	
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
        return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        secondOk = false;
        is_list = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String in1 = m_infile1.getStringValue();
    	String in2 = m_infile2.getStringValue();
    	
    	secondOk = false;
    	is_list = false;
    	
    	int end = checkEnding(in1);
    	//check first input file
    	if(Files.notExists(Paths.get(in1)) || in1.equals("")) {
    		throw new InvalidSettingsException("First input file does not exist");
    	} else if(end==0) {
    		this.setWarningMessage("The input file is not in a supported NGS format");
    	} else if(end==-1) {
    		is_list = true;
			String firstLine;
			try {
				firstLine = IO.head(Paths.get(in1), 1).get(0);
				in1 = firstLine.trim();
				String [] sep = {"\t", "\\s", ";", ","};
				for(String s:sep) {
					if(firstLine.contains(s)) {
						String [] col = firstLine.split(s);
						if(col.length > 2) setWarningMessage("Only the first or first two columns are read. Further columns are ignored.");
						in1 = col[0];
						in2 = col[1];
						this.sep = s;
						break;
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    	}
    	
    	//first input file is ok
    	dcs1 = new DataColumnSpecCreator(OUT_COL1, TYPES[checkEnding(in1)]).createSpec();
    	
    	if(in2.length()>0) {
    		String in2nozip = IO.removeZipExtension(in2);
    		if(Files.notExists(Paths.get(in2)) || !(in2nozip.endsWith(".fastq") || in2nozip.endsWith(".fq"))) {
    			setWarningMessage("Second input file does not exist or has disallowed ending and will be ignored!");
    		} else {
    			secondOk = true;
    			dcs2 = new DataColumnSpecCreator(OUT_COL2, TYPES[checkEnding(in2)]).createSpec();
    		}
    	}
    	
    	if(secondOk) {
    		specs = new DataColumnSpec[]{dcs1, dcs2};
    	} else {
    		specs = new DataColumnSpec[]{dcs1};
    	}
    	
    	return new DataTableSpec[]{new DataTableSpec(specs)};
    }
    
    /**
     * checks ending of input file
     * @param path 
     * @return -1 if input file is a list (.csv, .tsv, .list), the index of the ENDINGS or 0 if input file is not supported
     */
    protected int checkEnding(String path) {
    	if(path.endsWith(".csv") || path.endsWith(".tsv") || path.endsWith(".list")) {
    		return -1;
    	}
    	path = IO.removeZipExtension(path);
    	for(int i = 1; i < ENDINGS.length; i++) {
    		if(path.endsWith(ENDINGS[i])) {
    			return i;
    		}
    	}
    	return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_infile1.saveSettingsTo(settings);
         m_infile2.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_infile1.loadSettingsFrom(settings);
        m_infile2.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_infile1.validateSettings(settings);
        m_infile2.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}

