package de.helmholtz_muenchen.ibis.ngs.vcfSampler;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.BioEntity;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;
import de.helmholtz_muenchen.ibis.utils.ngs.VEPAnnotationParser;

/**
 * This is the model implementation of VCFSampler.
 * 
 *
 * @author Tim Jeske
 */
public class VCFSamplerNodeModel extends NodeModel {
    
	//configuration keys
	static final String CFGKEY_BUFFER = "buffer_size";
	static final String CFGKEY_CASES = "cases";
	static final String CFGKEY_CTRLS = "controls";
	static final String CFGKEY_DEF = "signal_def";
	static final String CFGKEY_NOISE = "noise_def";
	
	//settings models
	private final SettingsModelInteger m_buffer = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_BUFFER, 50);
	private final SettingsModelInteger m_cases = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CASES, 100);
	private final SettingsModelInteger m_ctrls = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CTRLS, 100);
	private final SettingsModelString m_def = new SettingsModelString(VCFSamplerNodeModel.CFGKEY_DEF,"10/2");
	private final SettingsModelString m_noise = new SettingsModelString(VCFSamplerNodeModel.CFGKEY_NOISE,"100/2");
	
	private int vcf_index;
	public static final String OUT_COL1 = "Path2SampledVCF";
	public static final String OUT_COL2 = "PedFile";
	
    /**
     * Constructor for the node model.
     */
    protected VCFSamplerNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	//input file
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(vcf_infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	
    	int buffer_size = m_buffer.getIntValue();
    	int case_size = m_cases.getIntValue();
    	int ctrls_size = m_ctrls.getIntValue();
    	int samples = case_size + ctrls_size;
    	
    	int iter = this.getAvailableFlowVariables().get("currentIteration").getIntValue();
    	String outfile = IO.replaceFileExtension(vcf_infile, "sampled_"+samples+"_Iter_"+iter+".vcf");
    	
    	ArrayList<double []> probs_list = new ArrayList<>();
    	ArrayList<String> var_fields = new ArrayList<>();
    	ArrayList<Set<String>> gene_list = new ArrayList<>();
    	VCFFile vcf_it = new VCFFile(vcf_infile);
    	AnnotationParser ap = new VEPAnnotationParser(vcf_it.getInfoHeader(VEPAnnotationParser.ANN_ID));
    	
    	writeHeader(outfile,vcf_it.getCompleteHeader());

    	//read input VCF
    	VCFVariant var;
//    	double ac_het, ac_hom, prob_het, prob_hom;
    	double an, ac, prob_un;
    	while(vcf_it.hasNext()) {
    		var = vcf_it.next();
    		
//    		ac_het = Double.parseDouble(var.getInfoField("AC_Het"));
//    		ac_hom = Double.parseDouble(var.getInfoField("AC_Hom"));
    		ac = Double.parseDouble(var.getInfoField("AC_Adj"));
    		an = Double.parseDouble(var.getInfoField("AN_Adj"));
    		
//    		prob_het = 2.0 * (ac_het/an);
//    		prob_hom = 2.0 * (ac_hom/an);
    		if(ac == 0.0 || var.getChrom().contains("Y") || var.getChrom().contains("X")) {
    			continue;
    		}
//    		prob_un = 1 - prob_het - prob_hom;
    		prob_un = Math.pow(1.0 - (ac/an),2);
//    		probs_list.add(new double[]{prob_hom+ prob_het,prob_un});
    		probs_list.add(new double[]{1-prob_un,prob_un});
    		var_fields.add(var.getChrom()+"\t"+var.getPos()+"\t"+var.getId()+"\t"+var.getRef()+"\t"+var.getAlt()+"\t"+var.getQual()+"\t"+var.getFilter()+"\t"+var.getInfo()+"\t"+var.getFormat());
    		gene_list.add(ap.getEntity2AlleleIds(var.getInfoField(ap.getAnnId()), BioEntity.GENE_ID).keySet());
    	}
    	
    	//create noise and document
    	String noise_def = m_noise.getStringValue();
    	String [] defs = noise_def.split(";");
    	int nr_rands, random;
    	double increase, bg;
    	double [] tmp;
    	
    	int range = probs_list.size();
    	HashSet<Integer> general_rand_indices = new HashSet<>();
    	HashSet<Integer> my_rand_indices;
    	String nl = System.getProperty("line.separator");
    	
    	StringBuilder sb = new StringBuilder();
    	
    	for(String d: defs) {
    		nr_rands = Integer.parseInt(d.split("/")[0]);
    		increase = Double.parseDouble(d.split("/")[1]);
    		
        	my_rand_indices = new HashSet<>();
        	while(my_rand_indices.size() < nr_rands) {
        		random = new Random().nextInt(range);
        		if(general_rand_indices.contains(random)) continue;
        		general_rand_indices.add(random);
        		my_rand_indices.add(random);
        	}
        	
        	for(int i : my_rand_indices) {
        		tmp = probs_list.get(i);
        		bg = tmp[0];
        		tmp[0] = increase*tmp[0];
        		if(tmp[0] > 1.0) {
        			tmp[0] = 1.0;
        		}
        		tmp[1] = 1 - tmp[0];
        		probs_list.set(i, tmp);
        		sb.append(gene_list.get(i)+"\t"+bg+"\t"+increase+"\t"+tmp[0]+nl);
        	}
    	}
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(IO.replaceFileExtension(outfile, ".noise.tsv")));
    	bw.write(sb.toString());
    	bw.close();
    	
    	
    	ArrayList<double []> case_list = new ArrayList<>();
    	double [] my;
    	for(double [] a: probs_list) {
    		my = new double[a.length];
    		for(int i = 0; i < a.length; i++) {
    			my[i] = a[i];
    		}
    		case_list.add(my);
    	}
    	
    	//create signals and document
    	String signal_def = m_def.getStringValue();
    	defs = signal_def.split(";");
    	sb = new StringBuilder();
    	
    	for(String d: defs) {
    		nr_rands = Integer.parseInt(d.split("/")[0]);
    		increase = Double.parseDouble(d.split("/")[1]);
    		
        	my_rand_indices = new HashSet<>();
        	while(my_rand_indices.size() < nr_rands) {
        		random = new Random().nextInt(range);
        		if(general_rand_indices.contains(random)) continue;
        		general_rand_indices.add(random);
        		my_rand_indices.add(random);
        	}
        	
        	for(int i : my_rand_indices) {
        		tmp = case_list.get(i);
        		bg = tmp[0];
        		tmp[0] = increase*tmp[0];
        		if(tmp[0] > 1.0) {
        			tmp[0] = 1.0;
        		}
        		tmp[1] = 1 - tmp[0];
        		case_list.set(i, tmp);
        		sb.append(gene_list.get(i)+"\t"+bg+"\t"+increase+"\t"+tmp[0]+nl);
        	}
    	}
    	
    	bw = Files.newBufferedWriter(Paths.get(IO.replaceFileExtension(outfile, ".signals.tsv")));
    	bw.write(sb.toString());
    	bw.close();
    	
    	//create vcf
    	int [] coding = {1,0};
    	Statistics stats = new Statistics();
    	
    	int end;
    	for(int i = 0; i < probs_list.size(); i+= buffer_size) {
    		end = i + buffer_size;
    		if(end >= probs_list.size()) {
    			end = probs_list.size();
    		}
    		int [][] gts_case = stats.getSamples(coding, case_list.subList(i, end), case_size);
    		int [][] gts_ctrl = stats.getSamples(coding, probs_list.subList(i,end), ctrls_size);
        	writeVCF(var_fields.subList(i, end),outfile,gts_case, gts_ctrl);
    	}
    	
    	stats.quit();
    	
    	String ped_file = IO.replaceFileExtension(outfile, ".ped");
    	writePEDFile(ped_file);
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile),
    			(FileCell) FileCellFactory.create(ped_file)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	
        return new BufferedDataTable[]{outTable};
    }
    
    private void writeVCF(List<String> vars, String outfile, int [][] gts_case, int [][] gts_ctrl) throws IOException {
    	
    	String [] gts_rep = {"0/0", "0/1"};
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile), StandardOpenOption.APPEND);
    	
    	for(int i = 0; i < gts_case.length; i++) {
    		bw.write(vars.get(i));
    		for(int j = 0; j < gts_case[i].length; j++) {
    			bw.write("\t"+gts_rep[gts_case[i][j]]);
    		}
    		for(int j = 0; j < gts_ctrl[i].length; j++) {
    			bw.write("\t"+gts_rep[gts_ctrl[i][j]]);
    		}
    		bw.newLine();
    	}
    	bw.close();
    }
    
    private void writeHeader(String outfile, String header) throws IOException {
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header.trim()+"\tFORMAT");
    	for(int i = 0; i < m_cases.getIntValue()+m_ctrls.getIntValue(); i++) {
    		bw.write("\tS"+i);
    	}
    	bw.newLine();
    	bw.close();
    }
    
    private void writePEDFile(String outfile) throws IOException {
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	
    	for(int i = 0; i < m_cases.getIntValue(); i++) {
    		bw.write("F"+i+"\tS"+i+"\t0\t0\t-1\t2");
    		bw.newLine();
    	} 
    	for(int i= m_cases.getIntValue(); i < m_ctrls.getIntValue()+m_cases.getIntValue(); i++) {
    		bw.write("F"+i+"\tS"+i+"\t0\t0\t-1\t1");
    		bw.newLine();
    	}
    	bw.close();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	vcf_index = -1;
    	for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
    	
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}

    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_ctrls.saveSettingsTo(settings);
         m_buffer.saveSettingsTo(settings);
         m_cases.saveSettingsTo(settings);
         m_def.saveSettingsTo(settings);
         m_noise.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_ctrls.loadSettingsFrom(settings);
    	m_buffer.loadSettingsFrom(settings);
    	m_cases.loadSettingsFrom(settings);
    	m_def.loadSettingsFrom(settings);
    	m_noise.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_ctrls.validateSettings(settings);
    	m_buffer.validateSettings(settings);
    	m_cases.validateSettings(settings);
    	m_def.validateSettings(settings);
    	m_noise.validateSettings(settings);
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

