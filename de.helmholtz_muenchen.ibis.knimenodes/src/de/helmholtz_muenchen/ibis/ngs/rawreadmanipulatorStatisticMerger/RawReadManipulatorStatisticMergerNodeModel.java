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
package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeModel;


/**
 * This is the model implementation of RawReadManipulatorStatisticMergerNode.
 * 
 *
 * @author Michael Kluge
 */
public class RawReadManipulatorStatisticMergerNodeModel extends StatisticMergerNodeModel {
	
	public static final String END_FILENAME = ".filtered.stdOut.log";
	public static final String DATA_FILE = ".*" + END_FILENAME;
	public static final String NEWLINE = System.lineSeparator();
	public static final String TAB = "\t";
	public static final String FILE_ENDING = ".txt";
	protected static final String NAME_OF_FILE = "Filename";
	protected static final String MERGER_NAME = "RawReadManipulator";
	protected static final String MODULE_START = "Stats for Module: >(.*)<";
	protected static final String MODULE_SEP = " -> " + TAB;
			
	HashMap<String, HashMap<String, String>> VALUES = new HashMap<String, HashMap<String, String>>();

	static {
		// add module name
		ArrayList<String> modules = new ArrayList<String>();
		modules.add("Module statistic");
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
			
			String line = null;
			@SuppressWarnings("unused")
			String nameOfModule = null;
			String valueName = null;
			String value = null;
			boolean nextIsValue = false;
			Pattern p = Pattern.compile(MODULE_START);
			Matcher m = null;
			String filename = f.getName().replaceFirst(END_FILENAME, "");
			
			// read lines
			while((line = r.readLine()) != null) {
				line = line.trim();
				
				// check, if line starts with a module
				if(nextIsValue == false) {
					m = p.matcher(line);
	
					if(m.matches()) {
						nextIsValue = true;
						nameOfModule = m.group(1);
					}
				}
				// extract value
				else {
					// try to split line
					String tmp[] = line.split(MODULE_SEP);
					if(tmp.length == 2) {
						valueName = tmp[0].replace(" ", "_");
						value = tmp[1];
						
						// check, if name is already stored in hashmap
						if(!this.VALUES.containsKey(valueName)) 
							this.VALUES.put(valueName, new HashMap<String, String>());
						
						// store the value
						this.VALUES.get(valueName).put(filename, value);
					}
					// reset module status
					nextIsValue = false;
				}
			}	
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
		// write header
		StringBuffer buff = new StringBuffer(NAME_OF_FILE + TAB);
		HashSet<String> filenames = new HashSet<String>();
		LinkedHashSet<String> names = new LinkedHashSet<String>();
		
		// get header and filenames
		for(Iterator<String> it = this.VALUES.keySet().iterator(); it.hasNext(); ) {
			String name = it.next();
			buff.append(name);
			if(it.hasNext())
				buff.append(TAB);
			
			// get all files
			filenames.addAll(this.VALUES.get(name).keySet());
			names.add(name);
		}
		buff.append(NEWLINE);
		
		// get for each file all vars
		for(Iterator<String> it = filenames.iterator(); it.hasNext(); ) {
			String file = it.next();
			buff.append(file + TAB);
			// get all values
			for(Iterator<String> itNames = names.iterator(); itNames.hasNext(); ) {
				String name = itNames.next();
				if(this.VALUES.get(name).containsKey(file))
					buff.append(this.VALUES.get(name).get(file));
				else
					buff.append("0");
				
				if(itNames.hasNext())
					buff.append(TAB);
			}
			
			// make newline if not last file
			if(it.hasNext())
				buff.append(NEWLINE);
		}
		// write stuff to file
		outfile.write(buff.toString());
		outfile.flush();
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