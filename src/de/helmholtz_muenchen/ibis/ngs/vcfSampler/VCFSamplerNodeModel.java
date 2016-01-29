package de.helmholtz_muenchen.ibis.ngs.vcfSampler;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

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

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;

/**
 * This is the model implementation of VCFSampler.
 * 
 *
 * @author Tim Jeske
 */
public class VCFSamplerNodeModel extends NodeModel {
    
	private int vcf_index;
	public static final String OUT_COL1 = "Path2SampledVCF";
	
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
    	
    	int sample_size = 100;
    	
    	String outfile = IO.replaceFileExtension(vcf_infile, "sampled.vcf");
    	
    	//create probs list
    	ArrayList<double []> probs_list = new ArrayList<>();
    	VCFFile vcf_it = new VCFFile(vcf_infile);
    	VCFVariant var;
    	String ac_het, ac_hom, an;
    	double prob_un, prob_het, prob_hom;
    	while(vcf_it.hasNext()) {
    		var = vcf_it.next();
    		ac_het = var.getInfoField("AC_Het");
    		ac_hom = var.getInfoField("AC_Hom");
    		an = var.getInfoField("AN");
    		prob_het = Double.parseDouble(ac_het)/Double.parseDouble(an);
    		prob_hom = Double.parseDouble(ac_hom)/Double.parseDouble(an);
    		prob_un = 1 - prob_het - prob_hom;
    		probs_list.add(new double[]{prob_un, prob_het,prob_hom});
    	}
    	
    	//create double[][] array
    	double[][] probs = new double[probs_list.size()][3];
    	for(int i = 0; i < probs_list.size(); i++) {
    		probs[i] = probs_list.get(i);
    	}
    	
    	//create vcf
    	String [] elements = {"0/0", "0/1", "1/1"};
    	Statistics stats = new Statistics();
    	
    	String [][] gts = stats.getSamples(elements, probs, sample_size);
    	stats.quit();
    	
    	writeVCF(outfile,gts);
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	
        return new BufferedDataTable[]{outTable};
    }
    
    private void writeVCF(String outfile, String [][] gts) throws IOException {
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	
    	for(int i = 0; i < gts.length; i++) {
    		bw.write(gts[i][0]);
    		for(int j = 1; j < gts[i].length; j++) {
    			bw.write("\t"+gts[i][j]);
    		}
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
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
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

