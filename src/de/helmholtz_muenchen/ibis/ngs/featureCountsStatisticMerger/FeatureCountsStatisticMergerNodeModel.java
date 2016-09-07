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
package de.helmholtz_muenchen.ibis.ngs.featureCountsStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeModel;


/**
 * This is the model implementation of FeatureCountsStatisticMergerNode.
 * 
 *
 * @author Michael Kluge
 */
public class FeatureCountsStatisticMergerNodeModel extends StatisticMergerNodeModel {
	
	public static final String END_FILENAME = ".featureCounts.summary";
	public static final String DATA_FILE = ".*" + END_FILENAME;
	public static final String NEWLINE = System.lineSeparator();
	public static final String TAB = "\t";
	public static final String FILE_ENDING = ".txt";
	protected static final String NAME_OF_FILE = "Filename";
	protected static final String MERGER_NAME = "FeatureCounts";
	
	private final StringBuffer CONTENT = new StringBuffer();
	
	static {
		// add module name
		ArrayList<String> modules = new ArrayList<String>();
		modules.add("Assignment statistic");
		StatisticMergerNodeModel.addModuleNames(MERGER_NAME, modules);
	}
	
	@Override
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().matches(DATA_FILE))
			return false;
		
		// open file
		try {
			BufferedReader r = new BufferedReader(new FileReader(f));
			
			StringBuffer headerB = new StringBuffer();
			StringBuffer valuesB = new StringBuffer();
			
			String line = null;
			String name = null;
			String value = null;
			
			// read lines
			while((line = r.readLine()) != null) {
				// try to split line
				String tmp[] = line.split(TAB);
				if(tmp.length == 2) {
					// replace white spaces at start and end
					name = tmp[0].trim();
					value = tmp[1].trim();
										
					// if header should be written
					if(writeHeader) {
						headerB.append(name);
						headerB.append(TAB);
					}
					valuesB.append(value);
					valuesB.append(TAB);
				}
			}
			// add the filename
			headerB.append(NAME_OF_FILE);
			valuesB.append(f.getName().replaceFirst(END_FILENAME, ""));	

			// write to buffer
			if(writeHeader) 
				CONTENT.append(headerB.toString());
			
			CONTENT.append(NEWLINE);
			CONTENT.append(valuesB.toString());
			
			// close the file
			r.close();
			return true;
		}
		catch(Exception e) {
			e.printStackTrace(); 
		}
		return false;
	}
	
	@Override
	public void finalize(BufferedWriter outfile) throws IOException {
		outfile.write(CONTENT.toString());
		outfile.flush();
		this.CONTENT.delete(0, this.CONTENT.length());
	}

	@Override
	public String getNameOfStatisticFile() {
		return DATA_FILE;
	}

	@Override
	public String getMergerName() {
		return MERGER_NAME;
	}
}