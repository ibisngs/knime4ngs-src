/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.segemehl;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * This is the model implementation of Segemehl.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiter
 */
public class SegemehlNodeModel extends HTExecutorNodeModel {
	
	public static final String CFGKEY_SEGEMEHLFILE = "segemehlfile";
	public static final String CFGKEY_ACCURACY= "segemehAccuracy";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_READTYPE = "readtype";
	public static final String CFGKEY_THREADS = "threads";
	public static final String CFGKEY_CLIP5ADAPTER = "clip5adapter";
	public static final String CFGKEY_CLIP3ADAPTER = "clip3adapter";
	public static final String CFGKEY_ADAPTER5SEQ = "adapter5seq";
	public static final String CFGKEY_ADAPTER3SEQ = "adapter3seq";
	public static final String CFGKEY_AUTOADAPTER3SEQ = "autoadapter3seq";
	public static final String CFGKEY_CLIPPOLYA = "clippolya";
	public static final String CFGKEY_CLIPPINGACCURACY = "clippingaccuracy";
	public static final String CFGKEY_SOFTHARDCLIPPING = "softhardclipping";
//	public static final String CFGKEY_CHECKINDEX = "checkIndexRefSeq";
	public static final String CFGKEY_CHECKSPLITREADMAPPING = "checkSplitReadMapping";
	public static final String CFGKEY_CHECKSBISULFITEMAPPING = "checkBisulfiteMapping";
	public static final String CFGKEY_BISULFITEMAPPINGTYPE = "bisulfiteMappingType";
	public static final String CFGKEY_OPTIONAL = "optional";

	private final SettingsModelString m_segemehlfile = new SettingsModelString(SegemehlNodeModel.CFGKEY_SEGEMEHLFILE,"");
	private final SettingsModelString m_refseqfile = new SettingsModelString(SegemehlNodeModel.CFGKEY_REFSEQFILE,"");
	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_THREADS, 4, 1, 250);
	private final SettingsModelBoolean m_clip5adapter = new SettingsModelBoolean(CFGKEY_CLIP5ADAPTER, false);
	private final SettingsModelBoolean m_clip3adapter = new SettingsModelBoolean(CFGKEY_CLIP3ADAPTER, false);
	private final SettingsModelBoolean m_autoadapter3seq = new SettingsModelBoolean(CFGKEY_AUTOADAPTER3SEQ, false);
	private final SettingsModelString m_adapter5seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER5SEQ,"");
	private final SettingsModelString m_adapter3seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER3SEQ,"");
	private final SettingsModelBoolean m_clippolya = new SettingsModelBoolean(CFGKEY_CLIPPOLYA, false);
	private final SettingsModelIntegerBounded m_clippingaccuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_CLIPPINGACCURACY, 70, 0, 100);
	private final SettingsModelString m_softhardclipping = new SettingsModelString(SegemehlNodeModel.CFGKEY_SOFTHARDCLIPPING,"");
//	private final SettingsModelBoolean m_checkIndexRefSeq = new SettingsModelBoolean(CFGKEY_CHECKINDEX, true);
	private final SettingsModelBoolean m_checkSplitReadMapping = new SettingsModelBoolean(CFGKEY_CHECKSPLITREADMAPPING, false);
	private final SettingsModelBoolean m_checkBisulfiteMapping = new SettingsModelBoolean(CFGKEY_CHECKSBISULFITEMAPPING, false);
	private final SettingsModelString m_bisulfiteMappingType = new SettingsModelString(SegemehlNodeModel.CFGKEY_BISULFITEMAPPINGTYPE,"");
	private final SettingsModelIntegerBounded m_accuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_ACCURACY, 90, 0, 100);
	private final SettingsModelOptionalString m_optional = new SettingsModelOptionalString(SegemehlNodeModel.CFGKEY_OPTIONAL, "",false);

	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SAMFile";
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(SegemehlNodeModel.class);
	
	private String segemehl_bin, ref_genome;
	
    /**
     * Constructor for the node model.
     */
    protected SegemehlNodeModel() {
    
        super(1, 1, 2);
        
    	addSetting(m_segemehlfile);
    	addSetting(m_refseqfile);
    	addSetting(m_adapter3seq);
    	addSetting(m_adapter5seq);
    	addSetting(m_autoadapter3seq);
    	addSetting(m_clip3adapter);
    	addSetting(m_clip5adapter);
    	addSetting(m_clippolya);
    	addSetting(m_softhardclipping);
    	addSetting(m_threads);
    	addSetting(m_clippingaccuracy);
    	addSetting(m_checkSplitReadMapping);
    	addSetting(m_checkBisulfiteMapping);
    	addSetting(m_bisulfiteMappingType);
    	addSetting(m_accuracy);
    	addSetting(m_optional);
        
    	addPrefPageSetting(m_segemehlfile, IBISKNIMENodesPlugin.SEGEMEHL);
    	addPrefPageSetting(m_refseqfile, IBISKNIMENodesPlugin.REF_GENOME);
    	
        m_autoadapter3seq.setEnabled(false);
        m_adapter3seq.setEnabled(false);
        m_adapter5seq.setEnabled(false);
        m_clippingaccuracy.setEnabled(false);
        m_softhardclipping.setEnabled(false);
        m_bisulfiteMappingType.setEnabled(false);
//        m_checkBisulfiteMapping.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	int colIndex1 = 0;
    	if(SegemehlNodeDialog.getUseMainInputColBool()){
    		colIndex1 = inData[0].getDataTableSpec().findColumnIndex(SegemehlNodeDialog.getMainInputCol1());
    	}
     	
    	String path2reads1 = inData[0].iterator().next().getCell(colIndex1).toString();
    	String path2readFile2 	= "";
    	 
    	if(readType.equals("paired-end")){
    		int colIndex2 = 1;
        	if(SegemehlNodeDialog.getUseMainInputColBool()){
        		colIndex2 = inData[0].getDataTableSpec().findColumnIndex(SegemehlNodeDialog.getMainInputCol2());
        	}
    		path2readFile2 = inData[0].iterator().next().getCell(colIndex2).toString();	
    	}
    	  	
    	String path2indexedRefSeq = ref_genome.substring(0,ref_genome.lastIndexOf(".")+1)+"idx";
    	String basePath = path2reads1.substring(0,path2reads1.lastIndexOf('/')+1);
    	String outBaseName1 = path2reads1.substring(path2reads1.lastIndexOf("/")+1,path2reads1.lastIndexOf("."));
    	String outBaseName = outBaseName1;
    	String outBaseName2 = outBaseName1;
    	if(path2readFile2.length() > 1) {
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2reads1.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}
    	String outName = basePath+outBaseName+"_map.sam";
    	String outNameUnmatchedReads = basePath+outBaseName+"_unmatchedReads.f";
    	int nrOfThreads = m_threads.getIntValue();
    	int accuracy = m_accuracy.getIntValue();

    	
    // Indexing reference sequence: segemehl -x chr1.idx -d chr1.fa
    	if(Files.notExists(Paths.get(path2indexedRefSeq))) {
    		LOGGER.info("Indexing reference sequence.");
    		command.add(segemehl_bin);
    		command.add("-x "+path2indexedRefSeq);
	    	command.add("-d "+ref_genome);
	     	/**
	     	 * Execute
	     	 */
	    	File lockFile = new File(path2indexedRefSeq + ".klock");
			super.executeCommand(new String[]{StringUtils.join(command, " ")}, path2indexedRefSeq, exec, lockFile);
	     	command = new ArrayList<String>();	//Clear Array
	    	
    	} else {
    		LOGGER.info("Indexing reference sequence SKIPPED.");
    	}
    	
    	
    // Match/ align reads: segemehl.x -i chr1.idx -d chr1.fa -q myreads.fa > mymap.sam
    	LOGGER.info("Match/ align reads to indexed reference sequence.");
    	
    	command.add(segemehl_bin);
    	command.add("-s");
    	
    	if(m_clip5adapter.getBooleanValue()) {
    		command.add("-P " + m_adapter5seq.getStringValue());
    	}
    	if(m_clip3adapter.getBooleanValue()) {
//    		if(m_autoadapter3seq.getBooleanValue()) {
//    			command.add("-Y");
//    		} else {
    			command.add("-Q " + m_adapter3seq.getStringValue());
    		}
//    	}
    	if(m_clippolya.getBooleanValue()) {
    		command.add("-T");
    	}
    	if(m_clip5adapter.getBooleanValue() || m_clip3adapter.getBooleanValue() || m_clippolya.getBooleanValue()) {
    		command.add("-R " + m_clippingaccuracy.getIntValue());
    		if(m_softhardclipping.getStringValue().equals("Hard")) {
    			command.add("-C");
    		}
    	}
    	if(readType.equals("paired-end")) {
    		command.add("-p " + path2readFile2);
    	}
    	if(m_checkSplitReadMapping.getBooleanValue()) {
    		command.add("-S");
    	}
    	if(m_checkBisulfiteMapping.getBooleanValue()) {
    		String bmType = m_bisulfiteMappingType.getStringValue();
    		if(bmType.equals("methylC-seq/Lister et al.")) {
    			command.add("-F 1");
    		} else if(bmType.equals("bs-seq/Cokus et al. protocol")) {
    			command.add("-F 2");
    		} else if(bmType.equals("PAR-CLIP with 4SU")) {
    			command.add("-F 3");
    		} else if(bmType.equals("PAR-CLIP with 6SG")) {
    			command.add("-F 4");
    		}
    	}

    	command.add("-A "+accuracy);
    	command.add("-t "+nrOfThreads);
    	command.add("-i "+path2indexedRefSeq);
    	command.add("-d "+ref_genome);
    	command.add("-q "+path2reads1);
    	command.add("-o "+outName);
    	command.add("-u "+outNameUnmatchedReads);
    	
    	if(m_optional.isActive()){
    		command.add(m_optional.getStringValue());
    	}
    	
     	/**
     	 * Execute
     	 */
    	File lockFile = new File(outName + ".klock");
	
		// execute the command
		super.executeCommand(new String[]{StringUtils.join(command, " ")}, outName, exec, lockFile);


     	command = new ArrayList<String>();	//Clear Array
     	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			 FileCellFactory.create(outName)};
    	
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
    	
    	segemehl_bin = IO.processFilePath(m_segemehlfile.getStringValue());
    	ref_genome = IO.processFilePath(m_refseqfile.getStringValue());
  	
    	super.conf(inSpecs);
    	
		if(CompatibilityChecker.inputFileNotOk(segemehl_bin, false)) {
			throw new InvalidSettingsException("Set path to samtools binary!");
		}
    	
		
		if(ref_genome.length() > 1) {
			if(!FileValidator.checkFastaFormat(ref_genome)){
	            throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
	    	}
		}
    		
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, SAMCell.TYPE).createSpec()})};
    }
}

