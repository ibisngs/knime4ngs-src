package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulatorStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
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
			
	private final StringBuffer CONTENT = new StringBuffer();
	
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
			
			StringBuffer headerB = new StringBuffer();
			StringBuffer valuesB = new StringBuffer();
			
			String line = null;
			@SuppressWarnings("unused")
			String nameOfModule = null;
			String valueName = null;
			String value = null;
			boolean nextIsValue = false;
			Pattern p = Pattern.compile(MODULE_START);
			Matcher m = null;
			
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
						valueName = tmp[0];
						value = tmp[1];
						
						// if header should be written
						if(writeHeader) {
							headerB.append(valueName.replace(" ", "_"));
							headerB.append(TAB);
						}
						valuesB.append(value);
						valuesB.append(TAB);
					}
					// reset module status
					nextIsValue = false;
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