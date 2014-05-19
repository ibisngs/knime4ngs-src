package de.helmholtz_muenchen.ibis.ngs.starStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger.StatisticMergerNodeModel;


/**
 * This is the model implementation of StarStatisticMergerNode.
 * 
 *
 * @author Michael Kluge
 */
public class StarStatisticMergerNodeModel extends StatisticMergerNodeModel {
	
	public static final String DATA_FILE = "Log.final.out";
	public static final String NEWLINE = System.lineSeparator();
	public static final String TAB = "\t";
	public static final String FILE_ENDING = ".txt";
	protected static final String NAME_OF_FILE = "Filename";
	protected static final String MERGER_NAME = "Star";
	
	private static final HashMap<String, String> REPLACE_MAP = new HashMap<String, String>();
	private final StringBuffer CONTENT = new StringBuffer();
	private final String REGEX_NAME_END = " |";
	private final String REGEX_VALUES_REPLACE = "%";
	
	static {
		REPLACE_MAP.put("%", "percent");
		REPLACE_MAP.put(" ", "_");
		REPLACE_MAP.put(":", "_");
		REPLACE_MAP.put("/", "_");
		REPLACE_MAP.put(")", "");
		REPLACE_MAP.put("(", "");
		
		// add module name
		ArrayList<String> modules = new ArrayList<String>();
		modules.add("Mapping statistic");
		StatisticMergerNodeModel.addModuleNames(MERGER_NAME, modules);
	}
	
	@Override
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().equals(DATA_FILE))
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
					name = name.replace(REGEX_NAME_END, "");
					
					// replace non-allowed chars in name
					for(String rep : REPLACE_MAP.keySet()) {
						name = name.replace(rep, REPLACE_MAP.get(rep));
					}
					value = value.replaceAll(REGEX_VALUES_REPLACE, "");
					
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
			valuesB.append(f.getParentFile().getName());	

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