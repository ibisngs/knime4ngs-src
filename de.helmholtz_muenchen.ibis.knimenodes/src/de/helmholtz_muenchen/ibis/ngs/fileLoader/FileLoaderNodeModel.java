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
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastACell;
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
public class FileLoaderNodeModel extends SettingsStorageNodeModel {
	
	public static final String CFGKEY_INFILE1 = "infile1";
	public static final String CFGKEY_INFILE2 = "infile2";
	public static final String CFGKEY_ISLIST = "is_list";
	
	private final SettingsModelString m_infile1 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE1,"");
	private final SettingsModelString m_infile2 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE2,"");
	private final SettingsModelBoolean m_isList = new SettingsModelBoolean(FileLoaderNodeModel.CFGKEY_ISLIST,false);
	
	public static final String OUT_COL1 = "Path2File1";
	public static final String OUT_COL2 = "Path2File2";
	
	private static final String [] ENDINGS = {"",".fa", ".fasta",".g.vcf",".vcf",".fastq",".fq",".bam",".sam"};
	private static final DataType [] TYPES = {FileCell.TYPE, FastACell.TYPE, FastACell.TYPE, GVCFCell.TYPE, VCFCell.TYPE, FastQCell.TYPE, FastQCell.TYPE, BAMCell.TYPE, SAMCell.TYPE};
		
	private boolean secondOk = false;
	private String sep;
	private DataType type;
	
	private DataColumnSpec dcs1 = null;
	private DataColumnSpec dcs2 = null;
	private DataColumnSpec [] specs = null;
	
    /**
     * Constructor for the node model.
     */
    protected FileLoaderNodeModel() {
    	super(0, 1);
    	addSetting(m_infile1);
    	addSetting(m_infile2);
    	addSetting(m_isList);
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
    	
    	String in1 = IO.processFilePath(m_infile1.getStringValue());
    	String in2 = IO.processFilePath(m_infile2.getStringValue());
    	
    	if(m_isList.getBooleanValue()) {
    		BufferedReader bw = Files.newBufferedReader(Paths.get(in1));
    		String line;
    		while((line=bw.readLine())!=null) {
    			if(line.trim().equals("")) continue;
    			if(secondOk) {
    				String [] col = line.split(this.sep);
    				String file1 = IO.getAbsolutePath(new File(in1).getParent(), col[0]);
    				String file2 = IO.getAbsolutePath(new File(in1).getParent(), col[1]);
    				
    				
    				if(CompatibilityChecker.inputFileNotOk(file1, type)) {
    					setWarningMessage("Some input files are invalid!");
    					continue;
    				}
    				
    				if(CompatibilityChecker.inputFileNotOk(file2, type)) {
    					setWarningMessage("Some input files are invalid!");
    					continue;
    				}
    				
    				in1_list.add(file1);
    				if(col.length<2) {
    					setWarningMessage("Second column contains less entries than first! Ignoring those lines!");
    					continue;
    				}
    				in2_list.add(file2);
    			} else {
    				String file = IO.getAbsolutePath(new File(in1).getParent(), line.trim());
    				if(CompatibilityChecker.inputFileNotOk(file, type)) {
    					setWarningMessage("Some input files are invalid!");
    					continue;
    				}
    				in1_list.add(file);
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
        			FileCellFactory.create(in1_list.get(i))};
        	
        	if(secondOk) {
        		fileCell = new FileCell[] {
            			FileCellFactory.create(in1_list.get(i)),
            			FileCellFactory.create(in2_list.get(i))};
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String in1 = IO.processFilePath(m_infile1.getStringValue());
    	String in2 = IO.processFilePath(m_infile2.getStringValue());
    	
    	secondOk = false;
    	int end = checkEnding(in1);
    	
    	//check first input file
    	if(CompatibilityChecker.inputFileNotOk(in1)) {
    		throw new InvalidSettingsException("First input file does not exist or is empty!");
    	}
    	
    	if(m_isList.getBooleanValue()) {
    		if(end > 0) setWarningMessage("The content of the input file is probably not a list of file paths!");
    		String firstLine;
			try {
				firstLine = IO.head(Paths.get(in1), 1).get(0);
				String [] sep = {"\t", " ", ";", ","};
				for(String s:sep) {
					if(firstLine.contains(s)) {
						String [] col = firstLine.split(s);
						if(col.length > 2) setWarningMessage("Only the first or first two columns are read. Further columns are ignored.");
						in2 = IO.getAbsolutePath(new File(in1).getParent(), col[1]);
						in1 = IO.getAbsolutePath(new File(in1).getParent(), col[0]);
						this.sep = s;
						break;
					} else {
						in1 = IO.getAbsolutePath(new File(in1).getParent(), firstLine.trim());
					}
				}
			} catch (IOException e) {
				throw new InvalidSettingsException(e.getMessage());
			}
    	} else if(end==0) {
    		this.setWarningMessage("The input file is not in a supported NGS format");
    	} 
    	
    	type = TYPES[checkEnding(in1)];
    	
    	//repeat check of first input file as it might be a path in a list
    	if(CompatibilityChecker.inputFileNotOk(in1, type)) {
    		throw new InvalidSettingsException("First input file "+in1+" in your list does not exist or is empty!");
    	}
    	
    	//first input file is ok
    	dcs1 = new DataColumnSpecCreator(OUT_COL1, TYPES[checkEnding(in1)]).createSpec();
    	
    	if(in2.length()>0) {
    		if(CompatibilityChecker.inputFileNotOk(in2, type) || !TYPES[checkEnding(in2)].toString().equals("FastQCell")) {
    			setWarningMessage("Second input file does not exist, has a disallowed ending or is invalid and will be ignored!");
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
    	
    	System.out.println(dcs1.getType().toString());
    	
    	return new DataTableSpec[]{new DataTableSpec(specs)};
    }
    
    /**
     * checks ending of input file
     * @param path 
     * @return -1 if input file is a list (.csv, .tsv, .list), the index of the ENDINGS or 0 if input file is not supported
     */
    protected int checkEnding(String path) {
//    	if(path.endsWith(".csv") || path.endsWith(".tsv") || path.endsWith(".list")) {
//    		return -1;
//    	}
    	path = IO.removeZipExtension(path);
    	for(int i = 1; i < ENDINGS.length; i++) {
    		if(path.endsWith(ENDINGS[i])) {
    			return i;
    		}
    	}
    	return 0;
    }
}

