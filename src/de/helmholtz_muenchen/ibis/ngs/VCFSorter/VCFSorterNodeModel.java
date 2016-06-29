package de.helmholtz_muenchen.ibis.ngs.VCFSorter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;

import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

/**
 * This is the model implementation of VCFSorter.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class VCFSorterNodeModel extends HTExecutorNodeModel {
    
	static final String CFGKEY_REFSEQFILE = "refseqfile";
	final SettingsModelString m_refseqfile = new SettingsModelString(VCFSorterNodeModel.CFGKEY_REFSEQFILE,"");
	

	//output col names
	
	public static final String OUT_COL1 = "Path2SortedVCF";
		
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected VCFSorterNodeModel() {
    
        super(1, 1);
        addSetting(m_refseqfile);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	/*
    	 * check input file
    	 */
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	
    	if(Files.notExists(Paths.get(vcf_infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String infile_warning = CheckUtils.checkSourceFile(vcf_infile);
    	if(infile_warning != null) {
    		setWarningMessage(infile_warning);
    	}
    	
    	String outfile = IO.replaceFileExtension(vcf_infile, ".sorted.vcf");
    	String lockfile = IO.replaceFileExtension(outfile,SuccessfulRunChecker.LOCK_ENDING);
    	
    	ArrayList<String> cmd = new ArrayList<>();
		cmd.add("perl");
		cmd.add(IO.getScriptPath()+"scripts/perl/vcfsorter.pl");
		cmd.add(IO.replaceFileExtension(m_refseqfile.getStringValue(), ".dict"));
		cmd.add(vcf_infile);
		cmd.add(outfile);
    	
		String[] cmd_array = new String[cmd.size()];
		for (int i = 0; i < cmd.size(); i++) {
			cmd_array[i] = cmd.get(i);
		}
    	
		
		super.executeCommand(cmd_array, exec, new File(lockfile));
		
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

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	if(CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell")!=0){
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}	
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }

}

