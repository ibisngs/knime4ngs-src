package de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.lofsummary.MatrixSummary;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of CaseControlAnalyzer.
 * 
 *
 * @author Tim Jeske
 */
public abstract class CaseControlAnalyzerNodeModel extends NodeModel {
    
	//model file configure keys
	static final String CFGKEY_MODEL_GENE_ID = "model_gene_id";
	static final String CFGKEY_FREQ = "freq";
	static final String CFGKEY_POP_SIZE = "pop_size";
	
	//summary file configure keys
//	static final String CFGKEY_SUMMARY_GENE_ID = "summary_gene_id";
//	static final String CFGKEY_CASE_COND = "case_cond";
//	static final String CFGKEY_CASE_NCOND = "case_ncond";
//	static final String CFGKEY_CONTROL_COND = "control_cond";
//	static final String CFGKEY_CONTROL_NCOND = "control_ncond";
//	static final String CFGKEY_IGNORE_UNOBSERVED = "ignore_unobserved";
	
	//model file settings models
    private final SettingsModelString m_model_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_MODEL_GENE_ID, "gene_id");
    private final SettingsModelString m_freq = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_FREQ, "variant_freq");
    private final SettingsModelIntegerBounded m_pop_size = new SettingsModelIntegerBounded(CaseControlAnalyzerNodeModel.CFGKEY_POP_SIZE, 1, 1, Integer.MAX_VALUE);
    
    //summary file settings models
//    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
//    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
//    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
//    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
//    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_NCOND, "un_ctrl");
//    private final SettingsModelBoolean m_ignore_unobserved = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_IGNORE_UNOBSERVED,false);
    
	protected static final NodeLogger logger = NodeLogger.getLogger(CaseControlAnalyzerNodeModel.class);
    
    protected CaseControlAnalyzerNodeModel() {
        super(OptionalPorts.createOPOs(2), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	/*based on gene_summary file*/
//    	String summary_file;
//    	HashMap<String, ContingencyTable> gene2table;
//    	summary_file = inData[0].iterator().next().getCell(0).toString();
//    	if(summary_file.contains("gene_summary") || summary_file.contains("gene_set_summary")) {
//    		gene2table = readSummaryFile(summary_file, m_summary_gene_id.getStringValue(), m_case_cond.getStringValue(), m_case_ncond.getStringValue(), m_control_cond.getStringValue(), m_control_ncond.getStringValue());
//    	} else {
//    		gene2table = readMatrixFile(summary_file);
//    	}
    	
    	/*first port = input matrix*/
    	String matrix_file = inData[0].iterator().next().getCell(0).toString();
    	MatrixSummary ms = new MatrixSummary(matrix_file);
    	
    	/*second port = background model file*/
    	HashMap<String, Double> gene2frequency;
    	String model_file = inData[1].iterator().next().getCell(0).toString();
    	gene2frequency = readModelFile(model_file, m_model_gene_id.getStringValue(), m_freq.getStringValue());
    	
    	performAnalysis(inData, exec, gene2frequency, m_pop_size.getIntValue(), ms);
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(getOutCol(), FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(getOutfile())};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	return new BufferedDataTable[]{outTable};
    }

	private HashMap<String, Double> readModelFile(String file, String gene_id_header, String freq_header) throws IOException, InvalidSettingsException {
    	HashMap<String, Double> result = new HashMap<>();
    	
    	int gene_index = -1;
    	int freq_index = -1;
    	
    	String line;
    	String [] fields;
    	BufferedReader br = Files.newBufferedReader(Paths.get(file));
    	
    	String header = br.readLine();
    	if(header.contains(gene_id_header) && header.contains(freq_header)) {
    		String [] cols = header.split("\t");
    		for(int i = 0; i < cols.length; i++) {
    			if(cols[i].equals(gene_id_header)) {
    				gene_index = i;
    			} else if (cols[i].equals(freq_header)) {
    				freq_index = i;
    			}
    		}
    	} else {
    		throw new InvalidSettingsException("Specified column headers are not found in model file!");
    	}
    	
    	while((line=br.readLine())!=null) {
    		fields = line.split("\t");
    		result.put(fields[gene_index], Double.parseDouble(fields[freq_index]));
    	}
    	br.close();
    	return result;
    }
    
//	private HashMap<String, ContingencyTable> readMatrixFile(String file) throws IOException {
//		HashMap<String, ContingencyTable> result = new HashMap<>();
//		
//		HashMap<String, HashSet<String>> gene2aff_cases = new HashMap<>();
//		HashMap<String, HashSet<String>> gene2un_cases = new HashMap<>();
//		HashMap<String, HashSet<String>> gene2aff_ctrls = new HashMap<>();
//		HashMap<String, HashSet<String>> gene2un_ctrls = new HashMap<>();
//		HashSet<String> cases = new HashSet<>();
//		HashSet<String> ctrls = new HashSet<>();
//		
//		BufferedReader br = Files.newBufferedReader(Paths.get(file));
//		
//		String [] header = br.readLine().split("\t");
//		for(int i = 1; i < header.length; i++) {
//			if(header[i].contains("case")) {
//				cases.add(header[i]);
//			} else if(header[i].contains("control")) {
//				ctrls.add(header[i]);
//			}
//		}
//		
//		String line;
//		String [] fields;
//		String gene;
//		HashSet<String> tmp;
//		while((line=br.readLine())!=null) {
//			fields = line.split("\t");
//			gene = fields[0].split("_")[1];
//			if(!gene2aff_ctrls.containsKey(gene)) {
//				gene2aff_ctrls.put(gene, new HashSet<>());
//			}
//			if(!gene2aff_cases.containsKey(gene)) {
//				gene2aff_cases.put(gene, new HashSet<>());
//			}
//			if(!gene2un_ctrls.containsKey(gene)) {
//				gene2un_ctrls.put(gene, new HashSet<>(ctrls));
//			}
//			if(!gene2un_cases.containsKey(gene)) {
//				gene2un_cases.put(gene, new HashSet<>(cases));
//			}
//			
//			for(int i = 1; i < fields.length; i++) {
//				if(fields[i].equals("het") || fields[i].equals("hom") || fields[i].equals("pot_comp")) {
//					if(header[i].contains("case")) {
//						tmp = gene2aff_cases.get(gene);
//						tmp.add(header[i]);
//						gene2aff_cases.put(gene, tmp);
//						gene2un_cases.get(gene).remove(header[i]);
//					} else if(header[i].contains("control")) {
//						tmp = gene2aff_ctrls.get(gene);
//						tmp.add(header[i]);
//						gene2aff_ctrls.put(gene, tmp);
//						gene2un_ctrls.get(gene).remove(header[i]);
//					}
//				} else if(fields[i].equals("undef")) {
//					if(header[i].contains("case")) {
//						gene2un_cases.get(gene).remove(header[i]);
//					} else if(header[i].contains("control")) {
//						gene2un_ctrls.get(gene).remove(header[i]);
//					}
//				}
//			}
//		}
//		br.close();
//		
//		int case_aff, control_aff;
//		for(String g: gene2un_ctrls.keySet()) {
//			case_aff = gene2aff_cases.get(g).size();
//			control_aff = gene2aff_ctrls.get(g).size();
//			if(m_ignore_unobserved.getBooleanValue() && case_aff == 0 && control_aff == 0) {
//    			continue;
//    		}
//			result.put(g, new ContingencyTable(case_aff,gene2un_cases.get(g).size(),control_aff,gene2un_ctrls.get(g).size()));
//		}
//		return result;
//	}
	
//    private HashMap<String, ContingencyTable> readSummaryFile(String file, String gene_id, String case_cond, String case_ncond, String control_cond, String control_ncond) throws IOException, InvalidSettingsException {
//    	HashMap<String, ContingencyTable> result = new HashMap<>();
//    	
//    	int gene_index = -1;
//    	int a_index = -1;
//    	int b_index = -1;
//    	int c_index = -1;
//    	int d_index = -1;
//    	
//    	String line;
//    	String [] fields;
//    	BufferedReader br = Files.newBufferedReader(Paths.get(file));
//    	
//    	String header = br.readLine();
//    	if(header.contains(gene_id) && header.contains(case_cond) && header.contains(case_ncond) && header.contains(control_cond) && header.contains(control_ncond)) {
//    		String [] cols = header.split("\t");
//    		for(int i = 0; i < cols.length; i++) {
//    			if(cols[i].equals(gene_id)) {
//    				gene_index = i;
//    			} else if(cols[i].equals(case_cond)) {
//    				a_index = i;
//    			} else if(cols[i].equals(case_ncond)) {
//    				b_index = i;
//    			} else if(cols[i].equals(control_cond)) {
//    				c_index = i;
//    			} else if(cols[i].equals(control_ncond)) {
//    				d_index = i;
//    			}
//    		}
//    	} else {
//    		throw new InvalidSettingsException("Specified column headers are not found in summary file!");
//    	}
//    	
//    	int case_aff, case_un, control_aff, control_un;
//    	
//    	while((line=br.readLine())!=null) {
//    		fields = line.split("\t");
//    		case_aff = Integer.parseInt(fields[a_index]);
//    		case_un = Integer.parseInt(fields[b_index]);
//    		control_aff = Integer.parseInt(fields[c_index]);
//    		control_un = Integer.parseInt(fields[d_index]);
//    		if(m_ignore_unobserved.getBooleanValue() && case_aff == 0 && control_aff == 0) {
//    			continue;
//    		}
//	    	result.put(fields[gene_index], new ContingencyTable(case_aff, case_un, control_aff, control_un));
//    	}
//    	
//    	br.close();
//    	return result;
//    }
    

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
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	checkSettingsModelString(m_model_gene_id);
        checkSettingsModelString(m_freq);
    	
//    	checkSettingsModelString(m_summary_gene_id);
//    	checkSettingsModelString(m_case_cond);
//    	checkSettingsModelString(m_case_ncond);
//    	checkSettingsModelString(m_control_cond);
//    	checkSettingsModelString(m_control_ncond);
    	
        return new DataTableSpec[]{new DataTableSpec(
        		new DataColumnSpec[]{
        				new DataColumnSpecCreator(getOutCol(), FileCell.TYPE).createSpec()})};
    }
    
    private void checkSettingsModelString (SettingsModelString m) throws InvalidSettingsException {
    	String tmp = m.getStringValue();
    	if(tmp.equals("") || tmp==null) {
    		throw new InvalidSettingsException("The "+ m.getKey()+" field must not be empty!");
    	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_model_gene_id.saveSettingsTo(settings);
         m_freq.saveSettingsTo(settings);
         m_pop_size.saveSettingsTo(settings);
//       m_summary_gene_id.saveSettingsTo(settings);
//       m_case_cond.saveSettingsTo(settings);
//       m_case_ncond.saveSettingsTo(settings);
//       m_control_cond.saveSettingsTo(settings);
//       m_control_ncond.saveSettingsTo(settings);
//       m_ignore_unobserved.saveSettingsTo(settings);
         saveExtraSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_model_gene_id.loadSettingsFrom(settings);
    	m_freq.loadSettingsFrom(settings);
    	m_pop_size.loadSettingsFrom(settings);
//    	m_summary_gene_id.loadSettingsFrom(settings);
//    	m_case_cond.loadSettingsFrom(settings);
//    	m_case_ncond.loadSettingsFrom(settings);
//    	m_control_cond.loadSettingsFrom(settings);
//    	m_control_ncond.loadSettingsFrom(settings);
//    	m_ignore_unobserved.loadSettingsFrom(settings);
    	loadExtraValidatedSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_model_gene_id.validateSettings(settings);
        m_freq.validateSettings(settings);
        m_pop_size.validateSettings(settings);
//      m_summary_gene_id.validateSettings(settings);
//      m_case_cond.validateSettings(settings);
//      m_case_ncond.validateSettings(settings);
//      m_control_cond.validateSettings(settings);
//      m_control_ncond.validateSettings(settings);
//      m_ignore_unobserved.validateSettings(settings);
        validateExtraSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    }
    
	protected abstract void performAnalysis(BufferedDataTable[] inData, ExecutionContext exec,
			HashMap<String, Double> gene2frequency, int pop_size, MatrixSummary ms) throws IOException, InvalidSettingsException;
	
	protected abstract String getOutCol();
	protected abstract String getOutfile();
	protected abstract void saveExtraSettingsTo(final NodeSettingsWO settings);
    protected abstract void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException;
    protected abstract void validateExtraSettings(final NodeSettingsRO settings) throws InvalidSettingsException;
}